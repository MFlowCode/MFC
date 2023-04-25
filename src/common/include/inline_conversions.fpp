#:def s_compute_speed_of_sound()
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c)

        !$acc routine seq

        real(kind(0d0)), intent(IN) :: pres
        real(kind(0d0)), intent(IN) :: rho, gamma, pi_inf
        real(kind(0d0)), intent(IN) :: H
        real(kind(0d0)), dimension(2), intent(IN) :: adv
        real(kind(0d0)), intent(IN) :: vel_sum
        real(kind(0d0)), intent(OUT) :: c

        real(kind(0d0)) :: blkmod1, blkmod2

        integer :: q
        
        c = 1.0
    end subroutine s_compute_speed_of_sound
#:enddef

