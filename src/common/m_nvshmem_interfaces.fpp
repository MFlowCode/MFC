!>
!! @file
!! @brief Contains module m_nvshmem_interfaces

!> @brief bind(C) interface declarations for the NVSHMEM device backend.
module m_nvshmem_interfaces

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    use, intrinsic :: iso_c_binding, only: c_double, c_float, c_int, c_long_long, c_ptr, c_size_t

    implicit none

    interface
        function mfc_nvshmem_halo_value_bytes() bind(C) result(bytes)

            import c_int
            integer(c_int) :: bytes

        end function mfc_nvshmem_halo_value_bytes

        function mfc_nvshmem_initialize(elements, fortran_comm, inbox, inbox_elements) bind(C) result(status)

            import c_int, c_ptr, c_size_t
            integer(c_size_t), value       :: elements
            integer(c_int), value          :: fortran_comm
            type(c_ptr), intent(out)       :: inbox
            integer(c_size_t), intent(out) :: inbox_elements
            integer(c_int)                 :: status

        end function mfc_nvshmem_initialize

        function mfc_nvshmem_exchange(send_buffer, elements, dst_pe, src_pe, send_channel, recv_channel, &
                                      & recv_offset) bind(C) result(status)
#ifdef MFC_SINGLE_PRECISION
            import c_float, c_int, c_size_t
            real(c_float), intent(in) :: send_buffer(*)
#else
            import c_double, c_int, c_size_t
            real(c_double), intent(in) :: send_buffer(*)
#endif
            integer(c_size_t), value       :: elements
            integer(c_int), value          :: dst_pe, src_pe, send_channel, recv_channel
            integer(c_size_t), intent(out) :: recv_offset
            integer(c_int)                 :: status

        end function mfc_nvshmem_exchange

        function mfc_nvshmem_complete_receive(src_pe, recv_channel) bind(C) result(status)

            import c_int
            integer(c_int), value :: src_pe, recv_channel
            integer(c_int)        :: status

        end function mfc_nvshmem_complete_receive

        ! Startup IB-patch neighborhood exchange with point-to-point puts over symmetric buffers, ordered by quiet and barrier.
        function mfc_nvshmem_ib_patch_counts_initialize(num_pes) bind(C) result(status)

            import c_int
            integer(c_int), value :: num_pes
            integer(c_int)        :: status

        end function mfc_nvshmem_ib_patch_counts_initialize

        function mfc_nvshmem_ib_patch_counts_publish(count) bind(C) result(status)

            import c_int
            integer(c_int), value :: count
            integer(c_int)        :: status

        end function mfc_nvshmem_ib_patch_counts_publish

        function mfc_nvshmem_ib_patch_count_get(dst, src_pe) bind(C) result(status)

            import c_int
            integer(c_int), intent(out) :: dst
            integer(c_int), value       :: src_pe
            integer(c_int)              :: status

        end function mfc_nvshmem_ib_patch_count_get

        function mfc_nvshmem_ib_patch_inbox_initialize(slot_bytes, slots) bind(C) result(status)

            import c_int, c_size_t
            integer(c_size_t), value :: slot_bytes
            integer(c_int), value    :: slots
            integer(c_int)           :: status

        end function mfc_nvshmem_ib_patch_inbox_initialize

        function mfc_nvshmem_ib_patch_put(src, bytes, dst_pe, slot) bind(C) result(status)

            import c_int, c_ptr, c_size_t
            type(c_ptr), value       :: src
            integer(c_size_t), value :: bytes
            integer(c_int), value    :: dst_pe, slot
            integer(c_int)           :: status

        end function mfc_nvshmem_ib_patch_put

        function mfc_nvshmem_ib_patch_inbox_get(dst, bytes, slot) bind(C) result(status)

            import c_int, c_ptr, c_size_t
            type(c_ptr), value       :: dst
            integer(c_size_t), value :: bytes
            integer(c_int), value    :: slot
            integer(c_int)           :: status

        end function mfc_nvshmem_ib_patch_inbox_get

        function mfc_nvshmem_ib_patch_sync() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_ib_patch_sync

        function mfc_nvshmem_retire() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_retire

#ifdef MFC_HALO_SINGLE_ROUND
        function mfc_nvshmem_initialize_26(elements, inbox, inbox_elements) bind(C) result(status)

            import c_int, c_size_t, c_ptr
            integer(c_size_t), value :: elements
            type(c_ptr)              :: inbox
            integer(c_size_t)        :: inbox_elements
            integer(c_int)           :: status

        end function mfc_nvshmem_initialize_26

        function mfc_nvshmem_configure_26(neighbor_count, dst_pe, src_points, dst_points, point_count, &
                                          & channel) bind(C) result(status)

            import c_int, c_long_long
            integer(c_int), value            :: neighbor_count
            integer(c_int), intent(in)       :: dst_pe(*), channel(*)
            integer(c_long_long), intent(in) :: src_points(*), dst_points(*), point_count(*)
            integer(c_int)                   :: status

        end function mfc_nvshmem_configure_26

        function mfc_nvshmem_exchange_26(send_buffer, v_size, recv_offset) bind(C) result(status)

            import c_int, c_size_t, c_long_long, c_float, c_double
#ifdef MFC_SINGLE_PRECISION
            real(c_float), intent(in) :: send_buffer(*)
#else
            real(c_double), intent(in) :: send_buffer(*)
#endif
            integer(c_long_long), value :: v_size
            integer(c_size_t)           :: recv_offset
            integer(c_int)              :: status

        end function mfc_nvshmem_exchange_26

        function mfc_nvshmem_complete_26() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_complete_26

        function mfc_nvshmem_retire_26() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_retire_26

        function mfc_nvshmem_exchange_26_post(send_buffer, v_size) bind(C) result(status)

            import c_int, c_long_long, c_float, c_double
#ifdef MFC_SINGLE_PRECISION
            real(c_float), intent(in) :: send_buffer(*)
#else
            real(c_double), intent(in) :: send_buffer(*)
#endif
            integer(c_long_long), value :: v_size
            integer(c_int)              :: status

        end function mfc_nvshmem_exchange_26_post

        function mfc_nvshmem_exchange_26_wait(recv_offset) bind(C) result(status)

            import c_int, c_size_t
            integer(c_size_t) :: recv_offset
            integer(c_int)    :: status

        end function mfc_nvshmem_exchange_26_wait

        !> Device-resident epoch exchange (MFC_HALO_EPOCH_DEVICE). No per-exchange value crosses the host, so the exchange is
        !! CUDA-graph capturable. <recv_base> binds the device address of the Fortran nvshmem_26_recv_base scalar that k_advance_26
        !! then writes on device.
        function mfc_nvshmem_bind_epoch_26(recv_base) bind(C) result(status)

            import c_int, c_long_long
            integer(c_long_long) :: recv_base
            integer(c_int)       :: status

        end function mfc_nvshmem_bind_epoch_26

        function mfc_nvshmem_exchange_26_dev(send_buffer, v_size) bind(C) result(status)

            import c_int, c_long_long, c_float, c_double
#ifdef MFC_SINGLE_PRECISION
            real(c_float), intent(in) :: send_buffer(*)
#else
            real(c_double), intent(in) :: send_buffer(*)
#endif
            integer(c_long_long), value :: v_size
            integer(c_int)              :: status

        end function mfc_nvshmem_exchange_26_dev

        function mfc_nvshmem_exchange_26_post_dev(send_buffer, v_size) bind(C) result(status)

            import c_int, c_long_long, c_float, c_double
#ifdef MFC_SINGLE_PRECISION
            real(c_float), intent(in) :: send_buffer(*)
#else
            real(c_double), intent(in) :: send_buffer(*)
#endif
            integer(c_long_long), value :: v_size
            integer(c_int)              :: status

        end function mfc_nvshmem_exchange_26_post_dev

        function mfc_nvshmem_exchange_26_wait_dev() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_exchange_26_wait_dev

        function mfc_nvshmem_complete_26_dev() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_complete_26_dev
#endif
#ifdef MFC_SIMULATION
        function mfc_nvshmem_initialize_igr(elements, inbox, inbox_elements) bind(C) result(status)

            import c_int, c_size_t, c_ptr
            integer(c_size_t), value :: elements
            type(c_ptr)              :: inbox
            integer(c_size_t)        :: inbox_elements
            integer(c_int)           :: status

        end function mfc_nvshmem_initialize_igr

        function mfc_nvshmem_configure_igr(active_count, dst_pe, face, src_off, dst_off, count) bind(C) result(status)

            import c_int, c_long_long
            integer(c_int), value            :: active_count
            integer(c_int), intent(in)       :: dst_pe(*), face(*)
            integer(c_long_long), intent(in) :: src_off(*), dst_off(*), count(*)
            integer(c_int)                   :: status

        end function mfc_nvshmem_configure_igr

        function mfc_nvshmem_exchange_igr(send_buffer, recv_offset) bind(C) result(status)

            import c_int, c_size_t, c_float, c_double
#ifdef MFC_SINGLE_PRECISION
            real(c_float), intent(in) :: send_buffer(*)
#else
            real(c_double), intent(in) :: send_buffer(*)
#endif
            integer(c_size_t) :: recv_offset
            integer(c_int)    :: status

        end function mfc_nvshmem_exchange_igr

        function mfc_nvshmem_complete_igr() bind(C) result(status)

            import c_int
            integer(c_int) :: status

        end function mfc_nvshmem_complete_igr
#endif

        subroutine mfc_nvshmem_finalize() bind(C)

        end subroutine mfc_nvshmem_finalize
    end interface
#endif
end module m_nvshmem_interfaces
