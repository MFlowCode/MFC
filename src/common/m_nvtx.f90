!>
!! @file
!! @brief Contains module m_nvtx

!> @brief NVIDIA NVTX profiling API bindings for GPU performance instrumentation
module m_nvtx

    use iso_c_binding

    implicit none

    integer, private :: col(7) = [int(Z'0000ff00'), int(Z'000000ff'), int(Z'00ffff00'), int(Z'00ff00ff'), int(Z'0000ffff'), &
                            & int(Z'00ff0000'), int(Z'00ffffff')]

    character(len=256), private :: tempName

    !> Per-rank wall-clock accumulation of every named range, enabled by phase_timing_wrt
    integer, parameter           :: nvtx_name_len = 64, nvtx_max_timers = 256
    logical                      :: nvtx_timing = .false.
    integer                      :: nvtx_num_timers = 0
    character(len=nvtx_name_len) :: nvtx_timer_names(nvtx_max_timers)
    integer(c_int64_t)           :: nvtx_timer_ticks(nvtx_max_timers) = 0, nvtx_timer_calls(nvtx_max_timers) = 0
    integer, private             :: depth = 0, stack_id(64)
    integer(c_int64_t), private  :: stack_tick(64)

    type, bind(C) :: nvtxEventAttributes
        integer(c_int16_t) :: version = 1
        integer(c_int16_t) :: size = 48  !
        integer(c_int)     :: category = 0
        integer(c_int)     :: colorType = 1    !< NVTX_COLOR_ARGB = 1
        integer(c_int)     :: color
        integer(c_int)     :: payloadType = 0  !< NVTX_PAYLOAD_UNKNOWN = 0
        integer(c_int)     :: reserved0
        integer(c_int64_t) :: payload          !< union uint,int,double
        integer(c_int)     :: messageType = 1  !< NVTX_MESSAGE_TYPE_ASCII = 1
        type(c_ptr)        :: message          !< ascii char
    end type nvtxEventAttributes

#if defined(MFC_GPU) && defined(__PGI)
    interface nvtxRangePush
        ! push range with custom label and standard color
        subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')

            use iso_c_binding

            character(kind=c_char, len=*), intent(in) :: name

        end subroutine nvtxRangePushA
        ! push range with custom label and custom color
        subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')

            use iso_c_binding

            import :: nvtxEventAttributes
            type(nvtxEventAttributes), intent(in) :: event

        end subroutine nvtxRangePushEx
    end interface nvtxRangePush

    interface nvtxRangePop
        subroutine nvtxRangePop() bind(C, name='nvtxRangePop')

        end subroutine nvtxRangePop
    end interface nvtxRangePop
#endif

contains

    !> Push a named NVTX range for GPU profiling, optionally with a color based on the given identifier.
    subroutine nvtxStartRange(name, id)

        character(kind=c_char, len=*), intent(in) :: name
        integer, intent(in), optional             :: id
        type(nvtxEventAttributes)                 :: event

        if (nvtx_timing) call s_push_timer(name)

#if defined(MFC_GPU) && defined(__PGI)
        tempName = trim(name) // c_null_char

        if (.not. present(id)) then
            call nvtxRangePush(tempName)
        else
            event%color = col(mod(id, 7) + 1)
            event%message = c_loc(tempName)
            call nvtxRangePushEx(event)
        end if
#endif

    end subroutine nvtxStartRange

    !> Pop the current NVTX range to end the GPU profiling region.
    subroutine nvtxEndRange

#if defined(MFC_GPU) && defined(__PGI)
        call nvtxRangePop
#endif

        if (nvtx_timing .and. depth > 0) call s_pop_timer()

    end subroutine nvtxEndRange

    !> Index of name in list, or 0 if absent
    pure integer function f_nvtx_find(list, name) result(idx)

        character(len=*), intent(in) :: list(:), name

        do idx = 1, size(list)
            if (list(idx) == name) return
        end do
        idx = 0

    end function f_nvtx_find

    !> Open a timed range, registering its name on first use (id 0 = untracked)
    subroutine s_push_timer(name)

        character(len=*), intent(in) :: name
        integer                      :: id

        id = f_nvtx_find(nvtx_timer_names(1:nvtx_num_timers), name)
        if (id == 0 .and. nvtx_num_timers < nvtx_max_timers) then
            nvtx_num_timers = nvtx_num_timers + 1
            nvtx_timer_names(nvtx_num_timers) = name
            id = nvtx_num_timers
        end if
        depth = depth + 1
        stack_id(depth) = id
        call system_clock(stack_tick(depth))

    end subroutine s_push_timer

    !> Close the innermost timed range and accumulate its elapsed ticks
    subroutine s_pop_timer()

        integer(c_int64_t) :: tick
        integer            :: id

        call system_clock(tick)
        id = stack_id(depth)
        if (id > 0) then
            nvtx_timer_ticks(id) = nvtx_timer_ticks(id) + tick - stack_tick(depth)
            nvtx_timer_calls(id) = nvtx_timer_calls(id) + 1
        end if
        depth = depth - 1

    end subroutine s_pop_timer

    !> Accumulated seconds of every registered range
    function f_nvtx_timer_seconds() result(sec)

        real(c_double)     :: sec(nvtx_num_timers)
        integer(c_int64_t) :: rate

        call system_clock(count_rate=rate)
        sec = real(nvtx_timer_ticks(1:nvtx_num_timers), c_double)/real(rate, c_double)

    end function f_nvtx_timer_seconds

end module m_nvtx
