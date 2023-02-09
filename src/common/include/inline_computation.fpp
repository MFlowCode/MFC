    !> Computes the bubble number density n from the conservative variables
    !! @param vftmp is the void fraction
    !! @param nRtmp is the bubble number  density times the bubble radii
    !! @param ntmp is the output number bubble density
#:def comp_n_from_cons()

    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp)
        !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
    
        real(kind(0.d0)) :: nR3
        integer :: i

        nR3 = 0d0
        do i = 1, nb
            nR3 = nR3 + weight(i)*(nRtmp(i)**3d0)
        end do

        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)

    end subroutine s_comp_n_from_cons

#:enddef comp_n_from_cons
