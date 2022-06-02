
module nvtx

  use iso_c_binding

  implicit none

  integer, private :: col(7) = [ &
    INT(Z'0000ff00'), INT(Z'000000ff'), INT(Z'00ffff00'), &
    INT(Z'00ff00ff'), INT(Z'0000ffff'), INT(Z'00ff0000'), &
    INT(Z'00ffffff') &
  ]

  character(len=256), private :: tempName

  type, bind(C) :: nvtxEventAttributes
    integer(C_INT16_T) :: version = 1
    integer(C_INT16_T) :: size = 48        !
    integer(C_INT)     :: category = 0
    integer(C_INT)     :: colorType = 1    ! NVTX_COLOR_ARGB = 1
    integer(C_INT)     :: color
    integer(C_INT)     :: payloadType = 0  ! NVTX_PAYLOAD_UNKNOWN = 0
    integer(C_INT)     :: reserved0
    integer(C_INT64_T) :: payload          ! union uint,int,double
    integer(C_INT)     :: messageType = 1  ! NVTX_MESSAGE_TYPE_ASCII = 1 
    type   (C_PTR)     :: message          ! ascii char
  end type

#ifdef _OPENACC

  interface nvtxRangePush
    ! push range with custom label and standard color
    subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
      use iso_c_binding

      character(kind=C_CHAR,len=*) :: name
    end subroutine

    ! push range with custom label and custom color
    subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
      use iso_c_binding

      import :: nvtxEventAttributes
      type(nvtxEventAttributes) :: event
    end subroutine
  end interface

  interface nvtxRangePop
    subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
    end subroutine
  end interface

  contains

#endif

  contains

  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional :: id
    type(nvtxEventAttributes):: event

#ifdef _OPENACC

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
      call nvtxRangePush(tempName)
    else
      event%color=col(mod(id,7)+1)
      event%message=c_loc(tempName)
      call nvtxRangePushEx(event)
    end if

#endif
  end subroutine

  subroutine nvtxEndRange
#ifdef _OPENACC
    call nvtxRangePop
#endif
  end subroutine

end module nvtx
