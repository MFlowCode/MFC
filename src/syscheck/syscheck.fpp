#:def LOG(*args)
#ifdef MFC_MPI
    if (rank == 0) then
#endif
        print *, ${','.join(args)}$
#ifdef MFC_MPI
    end if
#endif
#:enddef LOG

#:def MPIC(*args)
#ifdef MFC_MPI
    @:LOG("[TEST] MPI: ${','.join([ x.replace("'", '') for x in args ])}$")
    ${','.join([ x.replace("'", '') for x in args ])}$
    if (ierr /= MPI_SUCCESS) then
        print *, " -> Error: ", ierr
        stop ierr
    end if
#else
    @:LOG("[SKIP] MPI: ${','.join([ x.replace("'", '') for x in args ])}$")
#endif
#:enddef MPIC

#:def ACCC(*args)
#ifdef MFC_OpenACC
    @:LOG("[TEST] ACC: ${','.join([ x.replace("'", '') for x in args ])}$")
    ${','.join([ x.replace("'", '') for x in args ])}$
#else
    @:LOG("[SKIP] ACC: ${','.join([ x.replace("'", '') for x in args ])}$")
#endif
#:enddef ACCC

#:def MPI(*args)
#ifdef MFC_MPI
    ${','.join([ x.replace("'", '') for x in args ])}$
#endif
#:enddef MPI

#:def ACC(*args)
#ifdef MFC_OpenACC
    ${','.join([ x.replace("'", '') for x in args ])}$
#endif
#:enddef ACC

program syscheck

    @:MPI(use mpi)
    @:ACC(use openacc)

    implicit none

    integer :: ierr, rank = 0, nRanks = 1

    @:ACC(integer(acc_device_kind) :: devtype)
    @:ACC(integer :: i, num_devices)
    @:ACC(real(8), allocatable, dimension(:) :: arr)
    @:ACC(integer, parameter :: N = 100)

    @:MPIC(call mpi_init(ierr))
    @:MPIC(call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr))
    @:MPIC(call mpi_barrier(MPI_COMM_WORLD, ierr))
    @:MPIC(call assert(rank >= 0))
    @:MPIC(call mpi_comm_size(MPI_COMM_WORLD, nRanks, ierr))
    @:MPIC(call assert(nRanks > 0 .and. rank < nRanks))

    @:ACCC('devtype = acc_get_device_type()')
    @:ACCC('num_devices = acc_get_num_devices(devtype)')
    @:ACCC(call assert(num_devices > 0))
    @:ACCC(call acc_set_device_num(mod(rank, nRanks), devtype))
    @:ACCC(allocate(arr(1:N)))
    @:ACCC('!$acc enter data create(arr(1:N))')
    @:ACCC('!$acc parallel loop')
    @:ACC(do i = 1, N)
    @:ACC(arr(i) = i)
    @:ACC(end do)
    @:ACCC('!$acc update host(arr(1:N))')
    @:ACCC('!$acc exit data delete(arr)')

    @:MPIC(call mpi_barrier(MPI_COMM_WORLD, ierr))
    @:MPIC(call mpi_finalize(ierr))

    @:LOG("")
    @:LOG("Syscheck: PASSED.")

end program syscheck

subroutine assert(condition)

    use iso_fortran_env, only: output_unit, error_unit

    logical, intent(in) :: condition

    if (.not. condition) then
        call flush (int(output_unit))
        call flush (int(error_unit))
        stop 1
    end if

end subroutine assert
