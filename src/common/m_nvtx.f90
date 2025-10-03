module m_nvtx

    use iso_c_binding

    implicit none

    integer, private :: col(7) = [ &
                        int(Z'0000ff00'), int(Z'000000ff'), int(Z'00ffff00'), &
                        int(Z'00ff00ff'), int(Z'0000ffff'), int(Z'00ff0000'), &
                        int(Z'00ffffff') &
                        ]

    character(len=256), private :: tempName

    type, bind(C) :: nvtxEventAttributes
        integer(c_int16_t) :: version = 1
        integer(c_int16_t) :: size = 48        !
        integer(c_int) :: category = 0
        integer(c_int) :: colorType = 1    ! NVTX_COLOR_ARGB = 1
        integer(c_int) :: color
        integer(c_int) :: payloadType = 0  ! NVTX_PAYLOAD_UNKNOWN = 0
        integer(c_int) :: reserved0
        integer(c_int64_t) :: payload          ! union uint,int,double
        integer(c_int) :: messageType = 1  ! NVTX_MESSAGE_TYPE_ASCII = 1
        type(c_ptr) :: message          ! ascii char
    end type nvtxEventAttributes

#if defined(MFC_GPU) && defined(__PGI)

    interface nvtxRangePush
        ! push range with custom label and standard color
        subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
            use iso_c_binding

            character(kind=c_char, len=*), intent(IN) :: name
        end subroutine nvtxRangePushA

        ! push range with custom label and custom color
        subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
            use iso_c_binding

            import :: nvtxEventAttributes
            type(nvtxEventAttributes), intent(IN) :: event
        end subroutine nvtxRangePushEx
    end interface nvtxRangePush

    interface nvtxRangePop
        subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
        end subroutine nvtxRangePop
    end interface nvtxRangePop

#endif

contains

    subroutine nvtxStartRange(name, id)
        character(kind=c_char, len=*), intent(IN) :: name
        integer, intent(IN), optional :: id
        type(nvtxEventAttributes) :: event

#if defined(MFC_GPU) && defined(__PGI)

        tempName = trim(name)//c_null_char

        if (.not. present(id)) then
            call nvtxRangePush(tempName)
        else
            event%color = col(mod(id, 7) + 1)
            event%message = c_loc(tempName)
            call nvtxRangePushEx(event)
        end if

#endif
    end subroutine nvtxStartRange

    subroutine nvtxEndRange
#if defined(MFC_GPU) && defined(__PGI)
        call nvtxRangePop
#endif
    end subroutine nvtxEndRange

end module m_nvtx
