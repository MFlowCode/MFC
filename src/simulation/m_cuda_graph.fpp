!>
!! @file
!! @brief Contains module m_cuda_graph

!> @brief Replays communication-free RHS device work from a captured CUDA graph instead of relaunching each kernel. A caller
!! brackets the region with f_graph_action and f_graph_commit and supplies a signature covering values that vary between calls;
!! physics and launch topology fixed for the process are captured automatically. Inert unless graph capture is enabled.
#:include 'macros.fpp'

module m_cuda_graph

#ifdef MFC_CUDA_GRAPHS
    use iso_c_binding

    implicit none

    integer, parameter              :: graph_replay = 0   !< the graph ran; skip the region
    integer, parameter              :: graph_capture = 1  !< recording; run the region, then commit
    integer, parameter              :: graph_plain = 2    !< run the region as usual
    integer, parameter              :: graph_slot_rhs = 0
    integer, parameter              :: graph_slot_bubble_ee = 1
    integer, parameter              :: graph_slot_fused_kernel = 2
    integer(c_long_long), parameter :: graph_queue = ${ACC_QUEUE if ACC_QUEUE is not None else 0}$_c_long_long

    interface
        function mfc_graph_action(slot, signature, queue) bind(c, name='mfc_graph_action')

            import :: c_int, c_long_long
            integer(c_int), value       :: slot
            integer(c_long_long), value :: signature
            integer(c_long_long), value :: queue
            integer(c_int)              :: mfc_graph_action

        end function mfc_graph_action

        function mfc_graph_commit(slot, queue) bind(c, name='mfc_graph_commit')

            import :: c_int, c_long_long
            integer(c_int), value       :: slot
            integer(c_long_long), value :: queue
            integer(c_int)              :: mfc_graph_commit

        end function mfc_graph_commit

        subroutine mfc_graph_finalize() bind(c, name='mfc_graph_finalize')

        end subroutine mfc_graph_finalize
    end interface

contains

    integer function f_graph_action(slot, signature)

        integer, intent(in)         :: slot
        integer(kind=8), intent(in) :: signature

        f_graph_action = int(mfc_graph_action(int(slot, c_int), int(signature, c_long_long), graph_queue))

    end function f_graph_action

    integer function f_graph_commit(slot)

        integer, intent(in) :: slot

        f_graph_commit = int(mfc_graph_commit(int(slot, c_int), graph_queue))

    end function f_graph_commit

    subroutine s_graph_signature(signature, value)

        integer(kind=8), intent(inout) :: signature
        integer(kind=8), intent(in)    :: value

        signature = ieor(signature, value - 7046029254386353131_8 + shiftl(signature, 6) + shiftr(signature, 2))

    end subroutine s_graph_signature

    subroutine s_finalize_cuda_graph_module

        call mfc_graph_finalize()

    end subroutine s_finalize_cuda_graph_module
#endif
end module m_cuda_graph
