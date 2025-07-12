#:include 'macros.fpp'
!>
!! @file m_patch_helper.f90
!! @brief Contains module m_helper

module m_patch_helper

    use m_global_parameters

    private; 
    public :: f_cut_on, &
              f_cut_off

contains

    function f_cut_on(x, eps) result(fx)

        real(wp), intent(in) :: x, eps
        real(wp) :: fx

        fx = 1 - f_gx(x/eps)/(f_gx(x/eps) + f_gx(1 - x/eps))

    end function f_cut_on

    function f_cut_off(x, eps) result(fx)

        real(wp), intent(in) :: x, eps
        real(wp) :: fx

        fx = f_gx(x/eps)/(f_gx(x/eps) + f_gx(1 - x/eps))

    end function f_cut_off

    function f_gx(x) result(gx)

        real(wp), intent(in) :: x
        real(wp) :: gx

        if (x > 0) then
            gx = exp(-1._wp/x)
        else
            gx = 0._wp
        end if

    end function f_gx

end module m_patch_helper
