    !> Computes the bubble number density n from the conservative variables
    !! @param vftmp is the void fraction
    !! @param nRtmp is the bubble number  density times the bubble radii
    !! @param ntmp is the output number bubble density
#:def comp_n_from_cons(vftmp, nRtmp, ntmp, weight)
        !$acc routine seq
        ! real(kind(0.d0)), intent(IN) :: vftmp
        ! real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        ! real(kind(0.d0)), intent(OUT) :: ntmp
        ! real(kind(0.d0)), dimension(nb) :: weight
    
        ! real(kind(0.d0)) :: nR3
    
        call s_quad(${nRtmp}$**3.d0, nR3, ${weight}$)  !returns itself if NR0 = 1
        ${ntmp}$ = DSQRT((4.d0*pi/3.d0)*nR3/${vftmp}$)
#:enddef comp_n_from_cons
