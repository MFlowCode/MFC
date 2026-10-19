"""Generate NVSHMEM backends from developer-authored MPI/Fypp regions."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

try:
    from .nvshmem_ir import TransportRegion, parse_sources
except ImportError:
    from nvshmem_ir import TransportRegion, parse_sources


@dataclass(frozen=True)
class CommunicationShape:
    halo_channel_count: int
    halo_26_channel_count: int
    halo_26_max_neighbors: int
    igr_face_count: int
    with_halo_26: bool = False
    with_igr: bool = False


def _require_int(source: str, pattern: str, description: str) -> int:
    match = re.search(pattern, source, re.MULTILINE)
    if match is None:
        raise ValueError(f"Cannot derive {description} from developer MPI/Fypp source")
    return int(match.group(1))


def _require_operations(region: TransportRegion, required: set[str]):
    operations = {call.operation for call in region.calls}
    if not required <= operations:
        raise ValueError(f"{region.name} requires MPI operations {sorted(required)}, found {sorted(operations)}")


def derive_transport_manifest(paths: list[Path]) -> list[dict[str, object]]:
    """Derive optimized NVSHMEM transport shapes from annotated MPI/Fypp code."""
    regions = {region.name: region for region in parse_sources(paths)}
    for name, operations in {
        "halo_axis": {"SENDRECV"},
    }.items():
        if name not in regions:
            raise ValueError(f"Missing required developer MPI region: {name}")
        _require_operations(regions[name], operations)
    for name, operations in {
        "halo_26": {"IRECV", "ISEND", "WAITALL"},
        "igr": {"IRECV", "ISEND", "WAITALL"},
        # Classified but not lowered; the NVSHMEM translations live inline behind MFC_NVSHMEM
        # in their modules (m_start_up.fpp, m_ibm.fpp), MPI is the fallback.
        "ib_patch_exchange": {"IRECV", "ISEND", "WAITALL"},
        "ib_force_exchange": {"IRECV", "ISEND", "WAITALL"},
        "ib_patch_migration": {"IRECV", "ISEND", "WAITALL"},
    }.items():
        if name in regions:
            _require_operations(regions[name], operations)

    sources = {path.name: path.read_text() for path in paths}
    mpi_source = sources["m_mpi_common.fpp"]
    dimensions = _require_int(mpi_source, r"integer\s*::\s*buffer_counts\(1:(\d+)\),\s*buffer_count", "halo dimensions")
    have_halo_26 = "halo_26" in regions
    if have_halo_26:
        max_neighbors = _require_int(
            mpi_source,
            r"integer,\s*parameter\s*::\s*halo_26_max_neighbors\s*=\s*(\d+)",
            "26-neighbor capacity",
        )
        tag = re.search(
            r"tag\s*=\s*(\d+)\s*\+\s*\(dz\s*\+\s*1\)\s*\*\s*(\d+)\s*" r"\+\s*\(dy\s*\+\s*1\)\s*\*\s*(\d+)\s*\+\s*\(dx\s*\+\s*1\)",
            mpi_source,
        )
        if tag is None:
            raise ValueError("Cannot derive the 26-neighbor MPI tag mapping")
        tag_base, plane_stride, row_stride = (int(value) for value in tag.groups())
        if row_stride**2 != plane_stride:
            raise ValueError("Unsupported 26-neighbor MPI tag mapping")
        neighbor_channels = row_stride * plane_stride
        if neighbor_channels != max_neighbors + 1:
            raise ValueError("MPI neighbor capacity and tag space disagree")

    have_igr = "igr" in regions
    if have_igr:
        igr_last = _require_int(
            mpi_source,
            r"integer\(kind=8\),\s*private\s*::\s*igr_face_off\(0:(\d+)\),\s*igr_face_cnt\(0:\d+\)",
            "IGR face count",
        )
    call_signatures = {name: [[call.operation, *call.arguments] for call in region.calls] for name, region in regions.items() if region.mode == "device"}
    transports = [
        {
            "name": "halo",
            "topology": "linear",
            "groups": dimensions,
            "lanes": 2,
            "channels": 2 * dimensions,
            "max_peers": 1,
            "reciprocal": "explicit",
            "put": "nbi_on_stream",
            "slots": 2,
            "mpi_regions": {"halo_axis": call_signatures["halo_axis"]},
        },
    ]
    if have_halo_26:
        transports.append(
            {
                "name": "halo_26",
                "topology": "cartesian",
                "groups": dimensions,
                "lanes": row_stride,
                "channels": neighbor_channels,
                "max_peers": max_neighbors,
                "reciprocal": "reverse",
                "put": "nbi_on_stream",
                "slots": 2,
                "tag_base": tag_base,
                "mpi_regions": {"halo_26": call_signatures["halo_26"]},
            }
        )
    if have_igr:
        transports.append(
            {
                "name": "igr",
                "topology": "linear",
                "groups": dimensions,
                "lanes": 2,
                "channels": igr_last + 1,
                "max_peers": igr_last + 1,
                "reciprocal": "xor1",
                "put": "nbi_on_stream",
                "slots": 2,
                "mpi_regions": {"igr": call_signatures["igr"]},
            }
        )
    return transports


def _optional_call(region: TransportRegion, operation: str, buffer: str | None = None) -> tuple[str, ...] | None:
    matches = [call.arguments for call in region.calls if call.operation == operation and (buffer is None or call.arguments[0].split("(", 1)[0] == buffer)]
    unique = list(dict.fromkeys(matches))
    if len(unique) > 1:
        raise ValueError(f"Expected one {operation} signature in {region.name}, found {len(unique)}")
    return unique[0] if unique else None


def _call(region: TransportRegion, operation: str, buffer: str | None = None) -> tuple[str, ...]:
    signature = _optional_call(region, operation, buffer)
    if signature is None:
        raise ValueError(f"Expected one {operation} signature in {region.name}, found 0")
    return signature


def render_lowerings(paths: list[Path]) -> dict[str, str]:
    """Lower parsed MPI calls to the existing optimized NVSHMEM Fortran calls."""
    regions = {region.name: region for region in parse_sources(paths)}
    halo_wp = _call(regions["halo_axis"], "SENDRECV", "buff_send")
    manifest = {transport["name"]: transport for transport in derive_transport_manifest(paths)}

    halo_axis_exchange = f"""        #:call GPU_HOST_DATA(use_device_addr='[{halo_wp[0]}]', no_async_wait=True)
            nvshmem_status = mfc_nvshmem_exchange({halo_wp[0]}, int({halo_wp[1]}, c_size_t), {halo_wp[3]}, {halo_wp[8]}, &
                                                  & send_channel, recv_channel, nvshmem_recv_offset)
        #:endcall GPU_HOST_DATA
"""

    lowerings = {
        "nvshmem_halo_axis_exchange.fpp": f"""        send_channel = 2*(mpi_dir - 1) + send_tag
        recv_channel = 2*(mpi_dir - 1) + recv_tag
        call nvtxStartRange(\"RHS-COMM-SENDRECV-NVSHMEM\")
{halo_axis_exchange}        call nvtxEndRange
        if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM halo exchange failed.')
""",
        "nvshmem_halo_axis_complete.fpp": """        nvshmem_status = mfc_nvshmem_complete_receive(src_proc, recv_channel)
        if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM halo acknowledgement failed.')
""",
    }
    if "halo_26" in regions:
        halo_26_send = _call(regions["halo_26"], "ISEND")
        halo_26_buffer = halo_26_send[0].split("(", 1)[0]
        lowerings["nvshmem_halo_26_exchange.fpp"] = f"""        call nvtxStartRange(\"RHS-COMM-26-SENDRECV-NVSHMEM\")
        if (halo_26_epoch_device) then
            ! Device-resident epoch exchange (MFC_HALO_EPOCH_DEVICE). The epoch, its parity,
            ! and the receive base live in device memory, so the launches below carry no
            ! per-exchange value and a captured CUDA graph replays them unchanged. The lazy
            ! bind below runs at the first exchange, before the graph backend ever captures.
            if (.not. halo_26_epoch_bound) then
                halo_26_epoch_bound = .true.
                #:call GPU_HOST_DATA(use_device_addr='[nvshmem_26_recv_base]')
                    nvshmem_status = mfc_nvshmem_bind_epoch_26(nvshmem_26_recv_base)
                #:endcall GPU_HOST_DATA
                if (nvshmem_status /= 0) call s_mpi_abort('Unable to bind the NVSHMEM device-epoch receive base.')
            end if
            if (do_post) then
                if (do_wait) then
                    #:call GPU_HOST_DATA(use_device_addr='[{halo_26_buffer}]', no_async_wait=True)
                        nvshmem_status = mfc_nvshmem_exchange_26_dev({halo_26_buffer}(0), int(v_size, c_long_long))
                    #:endcall GPU_HOST_DATA
                    if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor exchange failed.')
                else
                    #:call GPU_HOST_DATA(use_device_addr='[{halo_26_buffer}]', no_async_wait=True)
                        nvshmem_status = mfc_nvshmem_exchange_26_post_dev({halo_26_buffer}(0), int(v_size, c_long_long))
                    #:endcall GPU_HOST_DATA
                    if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor exchange post failed.')
                end if
            else if (do_wait) then
                nvshmem_status = mfc_nvshmem_exchange_26_wait_dev()
                if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor exchange wait failed.')
            end if
        else if (do_post) then
            if (do_wait) then
                #:call GPU_HOST_DATA(use_device_addr='[{halo_26_buffer}]', no_async_wait=True)
                    nvshmem_status = mfc_nvshmem_exchange_26({halo_26_buffer}(0), int(v_size, c_long_long), nvshmem_26_recv_offset)
                #:endcall GPU_HOST_DATA
                if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor exchange failed.')
                nvshmem_26_recv_base = int(nvshmem_26_recv_offset, kind=8)
                $:GPU_UPDATE(device='[nvshmem_26_recv_base]')
            else
                #:call GPU_HOST_DATA(use_device_addr='[{halo_26_buffer}]', no_async_wait=True)
                    nvshmem_status = mfc_nvshmem_exchange_26_post({halo_26_buffer}(0), int(v_size, c_long_long))
                #:endcall GPU_HOST_DATA
                if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor exchange post failed.')
            end if
        else if (do_wait) then
            nvshmem_status = mfc_nvshmem_exchange_26_wait(nvshmem_26_recv_offset)
            if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor exchange wait failed.')
            nvshmem_26_recv_base = int(nvshmem_26_recv_offset, kind=8)
            $:GPU_UPDATE(device='[nvshmem_26_recv_base]')
        end if
        call nvtxEndRange
"""
        lowerings["nvshmem_halo_26_complete.fpp"] = """        if (halo_26_epoch_device) then
            nvshmem_status = mfc_nvshmem_complete_26_dev()
        else
            nvshmem_status = mfc_nvshmem_complete_26()
        end if
        if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor completion failed.')
"""
        lowerings["nvshmem_halo_26_channel.fpp"] = f"""            nvshmem_26_channel(neighbor_idx) = &
                & int(f_26_neighbor_tag(dx, dy, dz) - {manifest['halo_26']['tag_base']}, c_int)
"""
    if "igr" in regions:
        igr_send = _call(regions["igr"], "ISEND")
        igr_buffer = igr_send[0].split("(", 1)[0]
        lowerings["nvshmem_igr_exchange.fpp"] = f"""        call nvtxStartRange(\"IGR-COMM-SENDRECV-NVSHMEM\")
        #:call GPU_HOST_DATA(use_device_addr='[{igr_buffer}]', no_async_wait=True)
            nvshmem_status = mfc_nvshmem_exchange_igr({igr_buffer}(0), nvshmem_igr_recv_offset)
        #:endcall GPU_HOST_DATA
        if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM IGR exchange failed.')
        nvshmem_igr_recv_base = int(nvshmem_igr_recv_offset, kind=8)
        $:GPU_UPDATE(device='[nvshmem_igr_recv_base]')
        call nvtxEndRange
"""
        lowerings["nvshmem_igr_complete.fpp"] = """        nvshmem_status = mfc_nvshmem_complete_igr()
        if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM IGR completion failed.')
"""
    return lowerings


def load_transport_shape(manifest: Path) -> CommunicationShape:
    """Load the derived transport registry and reject non-NBI modes."""
    entries = json.loads(manifest.read_text())
    transports = {entry["name"]: entry for entry in entries}
    if len(transports) != len(entries):
        raise ValueError("Duplicate NVSHMEM transport name")

    required = {"halo"}
    optional = {"halo_26", "igr"}
    missing = required - transports.keys()
    extra = transports.keys() - required - optional
    if missing or extra:
        raise ValueError(f"Unsupported NVSHMEM transport set (missing={sorted(missing)}, extra={sorted(extra)})")

    expected_layouts = {
        "halo": ("linear", "explicit"),
        "halo_26": ("cartesian", "reverse"),
        "igr": ("linear", "xor1"),
    }
    for name, transport in transports.items():
        topology, reciprocal = expected_layouts[name]
        if transport["topology"] != topology or transport["reciprocal"] != reciprocal:
            raise ValueError(f"Unsupported {name} NVSHMEM transport mapping")
        if transport["put"] != "nbi_on_stream" or transport["slots"] != 2:
            raise ValueError(f"{name} must use stream-ordered NBI puts with two flow-control slots")
        expected_channels = transport["groups"] * transport["lanes"] if topology == "linear" else transport["lanes"] ** transport["groups"]
        if transport["channels"] != expected_channels:
            raise ValueError(f"{name} channel count does not match its Fypp topology")

    return CommunicationShape(
        transports["halo"]["channels"],
        transports["halo_26"]["channels"] if "halo_26" in transports else 0,
        transports["halo_26"]["max_peers"] if "halo_26" in transports else 1,
        transports["igr"]["channels"] if "igr" in transports else 1,
        with_halo_26="halo_26" in transports,
        with_igr="igr" in transports,
    )


def _public_definitions(source: str) -> set[str]:
    return set(
        re.findall(
            r'extern\s+"C"\s+(?:int|void)\s+(mfc_nvshmem_[a-z0-9_]+)\s*\([^;]*?\)\s*\{',
            source,
            re.DOTALL,
        )
    )


def render(template: str, shape: CommunicationShape, source_hash: str) -> str:
    replacements = {
        "generated_header": "\n".join(
            [
                "// AUTO-GENERATED from developer MPI/Fypp communication regions by",
                f"// src/nvshmem/gen_nvshmem.py (sha256 {source_hash}). Do not edit.",
            ]
        ),
        "halo_channel_count": str(shape.halo_channel_count),
        "halo_26_channel_count": str(shape.halo_26_channel_count),
        "halo_26_max_neighbors": str(shape.halo_26_max_neighbors),
        "igr_face_count": str(shape.igr_face_count),
        "with_halo_26": "1" if shape.with_halo_26 else "0",
        "with_igr": "1" if shape.with_igr else "0",
    }
    for name, value in replacements.items():
        template = template.replace(f"@{name}@", value)
    unresolved = re.findall(r"@[a-z0-9_]+@", template)
    if unresolved:
        raise ValueError(f"Unresolved NVSHMEM template fields: {', '.join(sorted(set(unresolved)))}")
    return template


def generate(manifest: Path, host_template: Path, device_template: Path) -> tuple[str, str]:
    shape = load_transport_shape(manifest)
    source_hash = hashlib.sha256(manifest.read_bytes()).hexdigest()[:16]
    host = render(host_template.read_text(), shape, source_hash)
    device = render(device_template.read_text(), shape, source_hash)
    return host, device


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--mpi-source", type=Path, action="append")
    parser.add_argument("--manifest-output", type=Path)
    parser.add_argument("--lowerings-dir", type=Path)
    parser.add_argument("--host-template", type=Path, required=True)
    parser.add_argument("--device-template", type=Path, required=True)
    parser.add_argument("--host-output", type=Path, required=True)
    parser.add_argument("--device-output", type=Path, required=True)
    args = parser.parse_args()

    if args.mpi_source:
        if args.manifest_output is None:
            parser.error("--manifest-output is required with --mpi-source")
        manifest = derive_transport_manifest(args.mpi_source)
        args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
        args.manifest_output.write_text(json.dumps(manifest, indent=2) + "\n")
        if args.lowerings_dir is not None:
            args.lowerings_dir.mkdir(parents=True, exist_ok=True)
            for name, content in render_lowerings(args.mpi_source).items():
                wrapped = "#:include 'case.fpp'\n#:include 'macros.fpp'\n" + content
                (args.lowerings_dir / name).write_text(wrapped)
        manifest_path = args.manifest_output
    elif args.manifest is not None:
        manifest_path = args.manifest
    else:
        parser.error("either --manifest or --mpi-source is required")

    host, device = generate(manifest_path, args.host_template, args.device_template)
    args.host_output.parent.mkdir(parents=True, exist_ok=True)
    args.device_output.parent.mkdir(parents=True, exist_ok=True)
    args.host_output.write_text(host)
    args.device_output.write_text(device)


if __name__ == "__main__":
    main()
