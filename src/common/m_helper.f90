!>
!! @file m_helper.f90
!! @brief Contains module m_helper
module m_helper

    ! Dependencies =============================================================
    
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    ! ==========================================================================

    implicit none

    real(kind(0d0)) :: cart_y, cart_z
    real(kind(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    private; public :: s_compute_finite_difference_coefficients, &
        s_comp_n_from_prim, &
        s_comp_n_from_cons, &
        s_quad

contains

    !>  The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, s_cc, fd_coeff_s, buff_size, &
                                                                fd_number, fd_order, offset_s)

        integer :: lB, lE !< loop bounds
        integer, intent(IN) :: q
        integer, intent(IN) :: buff_size, fd_number, fd_order
        type(int_bounds_info), optional, intent(IN) :: offset_s
        real(kind(0d0)), allocatable, dimension(:, :), intent(INOUT) :: fd_coeff_s

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        integer :: i !< Generic loop iterator

        if (present(offset_s)) then
            lB = -offset_s%beg
            lE = q + offset_s%end
            allocate (fd_coeff_s(-fd_number:fd_number, lb:lE))
        else
            lB = 0
            lE = q
            allocate (fd_coeff_s(-fd_number:fd_number, lb:lE))
        endif

        ! Computing the 1st order finite-difference coefficients
        if (fd_order == 1) then
            do i = lB, lE
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order == 2) then
            do i = lB, lE
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = lB, lE
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

        deallocate(fd_coeff_s)

    end subroutine s_compute_finite_difference_coefficients ! --------------

    !> Computes the bubble number density n from the primitive variables
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp, weight)
        !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: Rtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: R3
        real(kind(0.d0)), dimension(nb) :: weight

        call s_quad(Rtmp**3.d0, R3, weight)

        ntmp = (3.d0/(4.d0*pi))*vftmp/R3

    end subroutine s_comp_n_from_prim

    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp, weight)
        !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp  
        real(kind(0.d0)) :: nR3
        real(kind(0.d0)), dimension(nb) :: weight

        call s_quad(nRtmp**3.d0, nR3, weight)

        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)

    end subroutine s_comp_n_from_cons

    !> Computes the quadrature for polydisperse bubble populations
        !! @param func is the bubble dynamic variables for each bin
        !! @param mom is the computed moment
    subroutine s_quad(func, mom, weight)
        !$acc routine seq
        real(kind(0.d0)), dimension(nb), intent(IN) :: func
        real(kind(0.d0)), intent(OUT) :: mom
        real(kind(0.d0)), dimension(nb) :: weight

        mom = dot_product(weight, func)

    end subroutine s_quad

end module m_helper
