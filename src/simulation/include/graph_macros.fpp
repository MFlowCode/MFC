#! One framework wrapper for a capture-safe region. CUDA records the kernels;
#! physics routines neither construct graph nodes nor manage capture/replay state.
#! The signature covers launch arguments that change between calls.
#:def GPU_CAPTURE_REGION(code, slot, signature, enabled)
#ifdef MFC_CUDA_GRAPHS
    block
        integer :: region_graph_action

        region_graph_action = graph_plain
        if (${enabled}$) region_graph_action = f_graph_action(${slot}$, ${signature}$)
        if (region_graph_action /= graph_replay) then
#endif
            $:code
#ifdef MFC_CUDA_GRAPHS
            if (region_graph_action == graph_capture) then
                if (f_graph_commit(${slot}$) /= 0) call s_mpi_abort('CUDA-graph region capture failed.')
            end if
        end if
    end block
#endif
#:enddef

#! Opening and closing halves of GPU_CAPTURE_REGION, for a region that must keep the indentation of the code it wraps
#:def GPU_CAPTURE_BEGIN(slot, signature, enabled)
#ifdef MFC_CUDA_GRAPHS
    block
        integer :: region_graph_action

        region_graph_action = graph_plain
        if (${enabled}$) region_graph_action = f_graph_action(${slot}$, ${signature}$)
        if (region_graph_action /= graph_replay) then
#endif
        #:enddef

        #:def GPU_CAPTURE_END(slot)
#ifdef MFC_CUDA_GRAPHS
            if (region_graph_action == graph_capture) then
                if (f_graph_commit(${slot}$) /= 0) call s_mpi_abort('CUDA-graph region capture failed.')
            end if
        end if
    end block
#endif
#:enddef

#! Skip the work a replayed graph already performed
#:def GPU_GRAPH_SKIP_BEGIN(action)
#ifdef MFC_CUDA_GRAPHS
    if (${action}$ /= graph_replay) then
#endif
    #:enddef

    #:def GPU_GRAPH_SKIP_END()
#ifdef MFC_CUDA_GRAPHS
    end if
#endif
#:enddef
