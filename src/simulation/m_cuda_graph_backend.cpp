#include <cuda_runtime.h>
#include <openacc.h>

#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>

namespace {
constexpr int max_slots = 4;
constexpr int max_entries = 16;
constexpr int warmup_entries = 8;

struct Entry {
    std::int64_t signature;
    cudaGraphExec_t executable;
};

struct Slot {
    Entry entries[max_entries];
    int count;
    int entries_seen;
    int capturing;
    bool disabled;
};

Slot slots[max_slots];
int captures_open;
bool slots_ready;

void ensure_ready() {
    if (slots_ready) {
        return;
    }
    slots_ready = true;
    for (int s = 0; s < max_slots; ++s) {
        slots[s].count = 0;
        slots[s].entries_seen = 0;
        slots[s].capturing = -1;
        slots[s].disabled = false;
        for (int e = 0; e < max_entries; ++e) {
            slots[s].entries[e].signature = 0;
            slots[s].entries[e].executable = nullptr;
        }
    }
}

cudaStream_t queue_stream(long long queue) {
    return static_cast<cudaStream_t>(acc_get_cuda_stream(queue));
}

// Active-box scan staging. The scan's per-row partial boxes land in a page-locked,
// device-mapped buffer, and an event recorded behind the scan marks when they are
// visible. The deferred apply waits on that event rather than on the whole queue.
cudaEvent_t scan_event;
int* scan_part = nullptr;

void disable(Slot& slot, const char* what, cudaError_t status) {
    slot.disabled = true;
    slot.capturing = -1;
    cudaGetLastError();
    std::fprintf(stderr, "MFC: CUDA-graph capture disabled at %s: %s\n", what, cudaGetErrorString(status));
}

}

extern "C" {
void mfc_scan_record(long long queue) {
    if (scan_event == nullptr) {
        if (cudaEventCreateWithFlags(&scan_event, cudaEventDisableTiming) != cudaSuccess) {
            cudaGetLastError();
            std::fprintf(stderr, "MFC: scan event unavailable; the deferred apply waits the queue\n");
            return;
        }
    }
    cudaEventRecord(scan_event, queue_stream(queue));
}

// Returns 0 once the staged scan's result is visible. Without a usable event it
// synchronizes the whole queue.
int mfc_scan_wait(long long queue) {
    if (scan_event != nullptr) {
        const cudaError_t status = cudaEventSynchronize(scan_event);
        if (status == cudaSuccess) return 0;
        cudaGetLastError();
    }
    return cudaStreamSynchronize(queue_stream(queue)) == cudaSuccess ? 0 : -1;
}

// Page-locked, device-mapped partial-box buffer of <lanes> ints, mapped into OpenACC so
// the scan loop addresses it as a present array.
int mfc_scan_part_init(long long lanes) {
    if (scan_part != nullptr || lanes <= 0) return 0;
    const std::size_t bytes = static_cast<std::size_t>(lanes)*sizeof(int);
    if (cudaHostAlloc(reinterpret_cast<void**>(&scan_part), bytes, cudaHostAllocMapped) != cudaSuccess) {
        cudaGetLastError();
        scan_part = nullptr;
        return -1;
    }
    void* device = nullptr;
    if (cudaHostGetDevicePointer(&device, scan_part, 0) != cudaSuccess || device == nullptr) {
        cudaGetLastError();
        cudaFreeHost(scan_part);
        scan_part = nullptr;
        return -1;
    }
    acc_map_data(scan_part, device, bytes);
    return 0;
}

void* mfc_scan_part_ptr() {
    return scan_part;
}

int mfc_graph_action(int slot_id, long long signature, long long queue) {
    ensure_ready();
    if (slot_id < 0 || slot_id >= max_slots) {
        return 2;
    }

    Slot& slot = slots[slot_id];
    if (slot.disabled || captures_open > 0) {
        return 2;
    }

    for (int e = 0; e < slot.count; ++e) {
        if (slot.entries[e].signature == signature) {
            cudaError_t status = cudaGraphLaunch(slot.entries[e].executable, queue_stream(queue));
            if (status != cudaSuccess) {
                disable(slot, "cudaGraphLaunch", status);
                return 2;
            }
            return 0;
        }
    }

    ++slot.entries_seen;
    if (slot.entries_seen <= warmup_entries || slot.count >= max_entries) {
        return 2;
    }

    cudaError_t status = cudaStreamBeginCapture(queue_stream(queue), cudaStreamCaptureModeRelaxed);
    if (status != cudaSuccess) {
        disable(slot, "cudaStreamBeginCapture", status);
        return 2;
    }
    slot.capturing = slot.count;
    slot.entries[slot.count].signature = signature;
    ++captures_open;
    return 1;
}
int mfc_graph_commit(int slot_id, long long queue) {
    ensure_ready();
    if (slot_id < 0 || slot_id >= max_slots) {
        return 0;
    }

    Slot& slot = slots[slot_id];
    if (slot.capturing < 0) {
        return 0;
    }
    const int entry = slot.capturing;
    slot.capturing = -1;
    --captures_open;

    cudaStream_t stream = queue_stream(queue);
    cudaGraph_t graph = nullptr;
    cudaError_t status = cudaStreamEndCapture(stream, &graph);
    if (status != cudaSuccess || graph == nullptr) {
        disable(slot, "cudaStreamEndCapture", status);
        return 1;
    }
    std::size_t nodes = 0;
    cudaGraphGetNodes(graph, nullptr, &nodes);
    if (nodes == 0) {
        cudaGraphDestroy(graph);
        slot.disabled = true;
        cudaGetLastError();
        std::fprintf(stderr, "MFC: CUDA-graph capture disabled: the region recorded no work\n");
        return 0;
    }

    status = cudaGraphInstantiate(&slot.entries[entry].executable, graph, 0);
    cudaGraphDestroy(graph);
    if (status != cudaSuccess) {
        disable(slot, "cudaGraphInstantiate", status);
        return 1;
    }

    status = cudaGraphLaunch(slot.entries[entry].executable, stream);
    if (status != cudaSuccess) {
        cudaGraphExecDestroy(slot.entries[entry].executable);
        slot.entries[entry].executable = nullptr;
        disable(slot, "cudaGraphLaunch", status);
        return 1;
    }

    slot.count = entry + 1;
    return 0;
}

void mfc_graph_finalize() {
    if (scan_event != nullptr) {
        cudaEventDestroy(scan_event);
        scan_event = nullptr;
    }
    if (scan_part != nullptr) {
        acc_unmap_data(scan_part);
        cudaFreeHost(scan_part);
        scan_part = nullptr;
    }
    if (!slots_ready) {
        return;
    }
    for (int s = 0; s < max_slots; ++s) {
        for (int e = 0; e < slots[s].count; ++e) {
            if (slots[s].entries[e].executable != nullptr) {
                cudaGraphExecDestroy(slots[s].entries[e].executable);
                slots[s].entries[e].executable = nullptr;
            }
        }
        slots[s].count = 0;
    }
}

}
