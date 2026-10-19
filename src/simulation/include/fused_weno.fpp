!> Shared reconstruction helpers for the native fused kernels; order 1 loads states directly.
#:set MUSCL = recon_type == 2
#:set WENO_RADIUS = muscl_order - 1 if MUSCL else weno_order//2
#:set TENO7 = teno and WENO_RADIUS == 3
#:set WENO_STENCILS = 5 if TENO7 else WENO_RADIUS + 1
#! Cross-product weight combination for plain JS and mapped WENO. QBMM, chemistry and
#! continuum damage use the per-stencil form of s_weno.
#:set CROSS_JS = not MUSCL and not (wenoz or teno) and not (bubbles_euler and qbmm) and not chemistry &
    & and not (hypoelasticity and cont_damage)
#:if CROSS_JS
    #:set XPROD = ['*'.join(f'weight_fac({j})' for j in range(WENO_STENCILS) if j != i) for i in range(WENO_STENCILS)]
    #:set DX = [f'd_coef({1 + i} + {WENO_STENCILS}*(normal - coef_lo))' for i in range(WENO_STENCILS)]
    #:set TSUM = ' + '.join(f'{XPROD[i]}*{DX[i]}' for i in range(WENO_STENCILS))
#:endif
#:if WENO_RADIUS > 0
    #:if not MUSCL
        #:set CONSTANT_STENCIL = ('vm2 == vp0 .and. vm1 == vp0 .and. vp1 == vp0 .and. vp2 == vp0' &
    & if WENO_RADIUS == 2 else ('vm1 == vp0 .and. vp1 == vp0' if WENO_RADIUS == 1 else &
    & 'vm3 == vp0 .and. vm2 == vp0 .and. vm1 == vp0 .and. vp1 == vp0 .and. vp2 == vp0 .and. vp3 == vp0'))
        #:if TENO7
            subroutine s_fused_weno_combine(poly, weight_fac, normal, d_coef, coef_lo, value)

            #:else
                subroutine s_fused_weno_combine(dvd, weight_fac, center, normal, poly_coef, d_coef, coef_lo, coef_n, mapped, value)

                #:endif

                $:GPU_ROUTINE(function_name='s_fused_weno_combine', parallelism='[seq]')

                #:if TENO7
                    real(wp), intent(in)       :: poly(0:4), weight_fac(0:4), d_coef(*)
                    integer, value, intent(in) :: normal, coef_lo
                    real(wp), intent(out)      :: value
                    real(wp)                   :: alpha(0:4), omega(0:4)
                #:else
                    real(wp), intent(in) :: dvd(-${WENO_RADIUS}$:${WENO_RADIUS - 1}$), weight_fac(0:${WENO_RADIUS}$), &
                         & poly_coef(*), d_coef(*)
                    real(wp), value, intent(in) :: center
                    integer, value, intent(in)  :: normal, coef_lo, coef_n
                    logical, value, intent(in)  :: mapped
                    real(wp), intent(out)       :: value
                    real(wp)                    :: poly(0:${WENO_RADIUS}$), alpha(0:${WENO_RADIUS}$), &
                         & omega(0:${WENO_RADIUS}$)${', rcp' if CROSS_JS else ''}$
                #:endif
                integer :: stencil

                #:if not TENO7
                    #:if WENO_RADIUS == 1
                        poly(0) = center + poly_coef(1 + normal - coef_lo)*dvd(0)
                        poly(1) = center + poly_coef(1 + normal - coef_lo + coef_n)*dvd(-1)
                    #:elif WENO_RADIUS == 2
                        poly(0) = center + poly_coef(1 + normal - coef_lo)*dvd(1) + poly_coef(1 + normal - coef_lo + 3*coef_n) &
                             & *dvd(0)
                        poly(1) = center + poly_coef(1 + normal - coef_lo + coef_n)*dvd(0) + poly_coef(1 + normal - coef_lo &
                             & + 4*coef_n)*dvd(-1)
                        poly(2) = center + poly_coef(1 + normal - coef_lo + 2*coef_n)*dvd(-1) + poly_coef(1 + normal - coef_lo &
                             & + 5*coef_n)*dvd(-2)
                    #:else
                        #:for ROW in range(4)
                            poly(${3 - ROW}$) = center + poly_coef(1 + normal - coef_lo + ${ROW}$*coef_n)*dvd(${2 - ROW}$) &
                                 & + poly_coef(1 + normal - coef_lo + ${ROW + 4}$*coef_n)*dvd(${1 - ROW}$) + poly_coef(1 + normal &
                                 & - coef_lo + ${ROW + 8}$*coef_n)*dvd(${-ROW}$)
                        #:endfor
                    #:endif
                #:endif

                #:if CROSS_JS
                    ! Single precision keeps the per-stencil form, whose weights stay in range.
#ifdef MFC_SINGLE_PRECISION
                #:endif
                do stencil = 0, ${WENO_STENCILS - 1}$
                    #:if wenoz or teno
                        ! WENO-Z and TENO supply a ready alpha factor.
                        alpha(stencil) = d_coef(1 + stencil + ${WENO_STENCILS}$*(normal - coef_lo))*weight_fac(stencil)
                    #:else
                        ! WENO-JS and mapped divide by the squared smoothness indicator, as in s_weno.
                        alpha(stencil) = d_coef(1 + stencil + ${WENO_STENCILS}$*(normal - coef_lo))/weight_fac(stencil)
                    #:endif
                end do
                #:if not teno
                    if (mapped) then
                        omega = alpha/sum(alpha)
                        do stencil = 0, ${WENO_STENCILS - 1}$
                            alpha(stencil) = (d_coef(1 + stencil + ${WENO_STENCILS}$*(normal - coef_lo))*(1._wp + d_coef(1 &
                                  & + stencil + ${WENO_STENCILS}$*(normal - coef_lo)) - 3._wp*omega(stencil)) + omega(stencil) &
                                  & **2._wp)*(omega(stencil)/(d_coef(1 + stencil + ${WENO_STENCILS}$*(normal - coef_lo))**2._wp &
                                  & + omega(stencil)*(1._wp - 2._wp*d_coef(1 + stencil + ${WENO_STENCILS}$*(normal - coef_lo)))))
                        end do
                    end if
                #:endif
                ! Per-stencil normalization, as in s_weno.
                omega = alpha/sum(alpha)
                value = ${' + '.join(f'omega({i})*poly({i})' for i in range(WENO_STENCILS))}$
                #:if CROSS_JS
#else
                    if (mapped) then
                        ! Mapped combination with the cross products t_i = d_i P_i in place of the
                        ! per-stencil divisions; the shared rcp factor cancels in the final ratio.
                        rcp = 1._wp/(${TSUM}$)
                        #:for I in range(WENO_STENCILS)
                            omega(${I}$) = ${XPROD[I]}$*${DX[I]}$*rcp
                        #:endfor
                        #:for I in range(WENO_STENCILS)
                            alpha(${I}$) = (${DX[I]}$*(1._wp + ${DX[I]}$ - 3._wp*omega(${I}$)) + omega(${I}$)**2._wp) &
                                  & *(${XPROD[I]}$*${DX[I]}$)/(${DX[I]}$**2._wp + omega(${I}$)*(1._wp - 2._wp*${DX[I]}$))
                        #:endfor
                        value = (${' + '.join(f'alpha({i})*poly({i})' for i in range(WENO_STENCILS))}$) &
                                 & /(${' + '.join(f'alpha({i})' for i in range(WENO_STENCILS))}$)
                    else
                        ! Plain JS as value = sum_i d_i P_i poly_i / sum_i d_i P_i with P_i = prod_{j /= i} w_j.
                        value = (${' + '.join(f'{XPROD[i]}*({DX[i]}*poly({i}))' for i in range(WENO_STENCILS))}$)/(${TSUM}$)
                    end if
#endif
                #:endif

            end subroutine s_fused_weno_combine
        #:endif

        #:if MUSCL
            subroutine s_fused_weno_state(q, q_index, normal_stride, epsilon_muscl, q_left, q_right)

            #:elif TENO7
                subroutine s_fused_weno_state(q, q_index, normal, normal_stride, d_l, d_r, coef_lo, epsilon_weno, teno_cutoff, &
                                              & q_left, q_right)
                #:else
                    subroutine s_fused_weno_state(q, q_index, normal, normal_stride, poly_l, poly_r, d_l, d_r, beta_coef, &
                                                  & coef_lo, coef_n, uniform, mapped, &
                                                  & epsilon_weno${', teno_cutoff' if teno else ''}$, q_left, q_right)
                    #:endif

                    $:GPU_ROUTINE(function_name='s_fused_weno_state', parallelism='[seq]')

                    real(stp), intent(in) :: q(*)
                    #:if MUSCL
                        integer, value, intent(in)  :: q_index, normal_stride
                        real(wp), value, intent(in) :: epsilon_muscl
                        real(wp), intent(out)       :: q_left, q_right
                    #:elif TENO7
                        real(wp), intent(in)        :: d_l(*), d_r(*)
                        integer, value, intent(in)  :: q_index, normal, normal_stride, coef_lo
                        real(wp), value, intent(in) :: epsilon_weno, teno_cutoff
                        real(wp), intent(out)       :: q_left, q_right
                    #:else
                        real(wp), intent(in)        :: poly_l(*), poly_r(*), d_l(*), d_r(*), beta_coef(*)
                        integer, value, intent(in)  :: q_index, normal, normal_stride, coef_lo, coef_n
                        logical, value, intent(in)  :: uniform, mapped
                        real(wp), value, intent(in) :: epsilon_weno${', teno_cutoff' if teno else ''}$
                        real(wp), intent(out)       :: q_left, q_right
                    #:endif
                    #:if MUSCL
                        real(wp) :: vm1, vp0, vp1, slopeL, slopeR, slope
                    #:elif TENO7
                        real(wp) :: vm3, vm2, vm1, vp0, vp1, vp2, vp3, beta(0:4), weight_fac(0:4), poly(0:4)
                    #:elif WENO_RADIUS == 3
                        real(wp) :: vm3, vm2, vm1, vp0, vp1, vp2, vp3, dvd(-3:2), beta(0:3), weight_fac(0:3)
                    #:elif WENO_RADIUS == 2
                        real(wp) :: vm2, vm1, vp0, vp1, vp2, dvd(-2:1), beta(0:2), weight_fac(0:2)
                    #:else
                        real(wp) :: vm1, vp0, vp1, dvd(-1:0), beta(0:1), weight_fac(0:1)
                    #:endif
                    #:if teno or (wenoz and not mapped_weno)
                        real(wp) :: tau
                    #:endif
                    #:if teno
                        real(wp) :: teno_sum
                    #:endif
                    #:if not MUSCL
                        integer :: stencil
                    #:endif
                    #:if hllc_mp_weno and WENO_RADIUS == 2 and not MUSCL
                        real(wp) :: curvature_lo, curvature_mid, curvature_hi, limited_lo, limited_hi
                        real(wp) :: four_a_b, four_b_a, curvature_sign, bound_ul, bound_md, bound_lc, bound_min, bound_max
                    #:endif

                    #:if WENO_RADIUS == 3
                        vm3 = real(q(q_index - 3*normal_stride), wp)
                    #:endif
                    #:if WENO_RADIUS >= 2
                        vm2 = real(q(q_index - 2*normal_stride), wp)
                    #:endif
                    vm1 = real(q(q_index - normal_stride), wp)
                    vp0 = real(q(q_index), wp)
                    vp1 = real(q(q_index + normal_stride), wp)
                    #:if WENO_RADIUS >= 2
                        vp2 = real(q(q_index + 2*normal_stride), wp)
                    #:endif
                    #:if WENO_RADIUS == 3
                        vp3 = real(q(q_index + 3*normal_stride), wp)
                    #:endif

                    #:if MUSCL
                        slopeL = vp1 - vp0
                        slopeR = vp0 - vm1
                        slope = 0._wp
                        #:if muscl_lim == 0
                            slope = 5e-1_wp*(slopeL + slopeR)
                        #:else
                            if (slopeL*slopeR > epsilon_muscl) then
                                #:if muscl_lim == 1
                                    slope = min(abs(slopeL), abs(slopeR))
                                #:elif muscl_lim == 2
                                    slope = min(2._wp*abs(slopeL), 2._wp*abs(slopeR))
                                    slope = min(slope, 5e-1_wp*(abs(slopeL) + abs(slopeR)))
                                #:elif muscl_lim == 3
                                    slope = ((slopeL + slopeR)*slopeL*slopeR)/(slopeL**2._wp + slopeR**2._wp)
                                #:elif muscl_lim == 4
                                    slope = 2._wp*slopeL*slopeR/(slopeL + slopeR)
                                #:elif muscl_lim == 5
                                    slope = -1._wp*min(-min(2._wp*abs(slopeL), abs(slopeR)), -min(abs(slopeL), 2._wp*abs(slopeR)))
                                #:endif
                            end if
                            #:if muscl_lim in (1, 2, 5)
                                if (slopeL < 0._wp) slope = -slope
                            #:endif
                        #:endif
                        q_left = vp0 - (5.e-1_wp*slope)
                        q_right = vp0 + (5.e-1_wp*slope)
                    #:else
                        #:if not teno
                            if (${CONSTANT_STENCIL}$) then
                                q_left = vp0
                                q_right = vp0
                                return
                            end if
                        #:endif

                        #:if TENO7
                            beta(0) = 13._wp/12._wp*(vm1 - 2._wp*vp0 + vp1)**2 + (vm1 - vp1)**2/4._wp + epsilon_weno
                            beta(1) = 13._wp/12._wp*(vp0 - 2._wp*vp1 + vp2)**2 + (3._wp*vp0 - 4._wp*vp1 + vp2)**2/4._wp &
                                 & + epsilon_weno
                            beta(2) = 13._wp/12._wp*(vm2 - 2._wp*vm1 + vp0)**2 + (vm2 - 4._wp*vm1 + 3._wp*vp0)**2/4._wp &
                                 & + epsilon_weno
                            beta(3) = (vp0*(2107._wp*vp0 - 9402._wp*vp1 + 7042._wp*vp2 - 1854._wp*vp3) + vp1*(11003._wp*vp1 &
                                 & - 17246._wp*vp2 + 4642._wp*vp3) + vp2*(7043._wp*vp2 - 3882._wp*vp3) + vp3*(547._wp*vp3)) &
                                 & /240._wp + epsilon_weno
                            beta(4) = (vm3*(547._wp*vm3 - 3882._wp*vm2 + 4642._wp*vm1 - 1854._wp*vp0) + vm2*(7043._wp*vm2 &
                                 & - 17246._wp*vm1 + 7042._wp*vp0) + vm1*(11003._wp*vm1 - 9402._wp*vp0) + vp0*(2107._wp*vp0)) &
                                 & /240._wp + epsilon_weno
                        #:else
                            #:if WENO_RADIUS == 3
                                dvd(2) = vp3 - vp2
                            #:endif
                            #:if WENO_RADIUS >= 2
                                dvd(1) = vp2 - vp1
                            #:endif
                            dvd(0) = vp1 - vp0
                            dvd(-1) = vp0 - vm1
                            #:if WENO_RADIUS >= 2
                                dvd(-2) = vm1 - vm2
                            #:endif
                            #:if WENO_RADIUS == 3
                                dvd(-3) = vm2 - vm3
                            #:endif

                            #:if WENO_RADIUS == 1
                                beta(0) = beta_coef(1 + normal - coef_lo)*dvd(0)*dvd(0) + epsilon_weno
                                beta(1) = beta_coef(1 + normal - coef_lo + coef_n)*dvd(-1)*dvd(-1) + epsilon_weno
                            #:elif WENO_RADIUS == 2
                                if (uniform) then
                                    beta(0) = 13._wp/12._wp*(dvd(1) - dvd(0))**2 + 0.25_wp*(dvd(1) - 3._wp*dvd(0))**2 + epsilon_weno
                                    beta(1) = 13._wp/12._wp*(dvd(0) - dvd(-1))**2 + 0.25_wp*(dvd(0) + dvd(-1))**2 + epsilon_weno
                                    beta(2) = 13._wp/12._wp*(dvd(-1) - dvd(-2))**2 + 0.25_wp*(3._wp*dvd(-1) - dvd(-2))**2 &
                                         & + epsilon_weno
                                else
                                    beta(0) = beta_coef(1 + normal - coef_lo)*dvd(1)*dvd(1) + beta_coef(1 + normal - coef_lo &
                                         & + 3*coef_n)*dvd(1)*dvd(0) + beta_coef(1 + normal - coef_lo + 6*coef_n)*dvd(0)*dvd(0) &
                                         & + epsilon_weno
                                    beta(1) = beta_coef(1 + normal - coef_lo + coef_n)*dvd(0)*dvd(0) + beta_coef(1 + normal &
                                         & - coef_lo + 4*coef_n)*dvd(0)*dvd(-1) + beta_coef(1 + normal - coef_lo + 7*coef_n) &
                                         & *dvd(-1)*dvd(-1) + epsilon_weno
                                    beta(2) = beta_coef(1 + normal - coef_lo + 2*coef_n)*dvd(-1)*dvd(-1) + beta_coef(1 + normal &
                                         & - coef_lo + 5*coef_n)*dvd(-1)*dvd(-2) + beta_coef(1 + normal - coef_lo + 8*coef_n) &
                                         & *dvd(-2)*dvd(-2) + epsilon_weno
                                end if
                            #:else
                                #:for ROW in range(4)
                                    beta(${3 - ROW}$) = beta_coef(1 + normal - coef_lo + ${ROW}$*coef_n)*dvd(${2 - ROW}$) &
                                         & *dvd(${2 - ROW}$) + beta_coef(1 + normal - coef_lo + ${ROW + 4}$*coef_n) &
                                         & *dvd(${2 - ROW}$)*dvd(${1 - ROW}$) + beta_coef(1 + normal - coef_lo &
                                         & + ${ROW + 8}$*coef_n)*dvd(${2 - ROW}$)*dvd(${-ROW}$) + beta_coef(1 + normal - coef_lo &
                                         & + ${ROW + 12}$*coef_n)*dvd(${1 - ROW}$)*dvd(${1 - ROW}$) + beta_coef(1 + normal &
                                         & - coef_lo + ${ROW + 16}$*coef_n)*dvd(${1 - ROW}$)*dvd(${-ROW}$) + beta_coef(1 + normal &
                                         & - coef_lo + ${ROW + 20}$*coef_n)*dvd(${-ROW}$)*dvd(${-ROW}$) + epsilon_weno
                                #:endfor
                            #:endif
                        #:endif

                        #:if teno
                            tau = abs(beta(${4 if TENO7 else WENO_RADIUS}$) - beta(${3 if TENO7 else 0}$))
                            do stencil = 0, ${WENO_STENCILS - 1}$
                                weight_fac(stencil) = 1._wp + tau/beta(stencil)
                                weight_fac(stencil) = (weight_fac(stencil)**3)**2
                            end do
                            teno_sum = sum(weight_fac)
                            do stencil = 0, ${WENO_STENCILS - 1}$
                                if (weight_fac(stencil)/teno_sum < teno_cutoff) then
                                    weight_fac(stencil) = 0._wp
                                else
                                    weight_fac(stencil) = 1._wp
                                end if
                            end do
                        #:elif wenoz and not mapped_weno
                            tau = abs(beta(${WENO_RADIUS}$) - beta(0))
                            ! WENO-Z uses this factor in place of the JS squared-smoothness divisor.
                            weight_fac = 1._wp
                            if (tau /= 0._wp) then
                                do stencil = 0, ${WENO_RADIUS}$
                                    #:if WENO_RADIUS == 3
                                        weight_fac(stencil) = 1._wp + (tau/beta(stencil))**${int(wenoz_q)}$
                                    #:else
                                        weight_fac(stencil) = 1._wp + tau/beta(stencil)
                                    #:endif
                                end do
                            end if
                        #:else
                            do stencil = 0, ${WENO_RADIUS}$
                                weight_fac(stencil) = beta(stencil)*beta(stencil)
                            end do
                        #:endif

                        #:if TENO7
                            #:for SIDE in ['l', 'r']
                                #:if SIDE == 'l'
                                    poly(0) = (2._wp*vm1 + 5._wp*vp0 - vp1)/6._wp
                                    poly(1) = (11._wp*vp0 - 7._wp*vp1 + 2._wp*vp2)/6._wp
                                    poly(2) = (-vm2 + 5._wp*vm1 + 2._wp*vp0)/6._wp
                                    poly(3) = (25._wp*vp0 - 23._wp*vp1 + 13._wp*vp2 - 3._wp*vp3)/12._wp
                                    poly(4) = (vm3 - 5._wp*vm2 + 13._wp*vm1 + 3._wp*vp0)/12._wp
                                #:else
                                    poly(0) = (-vm1 + 5._wp*vp0 + 2._wp*vp1)/6._wp
                                    poly(1) = (2._wp*vp0 + 5._wp*vp1 - vp2)/6._wp
                                    poly(2) = (2._wp*vm2 - 7._wp*vm1 + 11._wp*vp0)/6._wp
                                    poly(3) = (3._wp*vp0 + 13._wp*vp1 - 5._wp*vp2 + vp3)/12._wp
                                    poly(4) = (-3._wp*vm3 + 13._wp*vm2 - 23._wp*vm1 + 25._wp*vp0)/12._wp
                                #:endif
                                call s_fused_weno_combine(poly, weight_fac, normal, d_${SIDE}$, coef_lo, &
                                                          & q_${'left' if SIDE == 'l' else 'right'}$)
                            #:endfor
                        #:else
                            call s_fused_weno_combine(dvd, weight_fac, vp0, normal, poly_l, d_l, coef_lo, coef_n, mapped, q_left)
                            call s_fused_weno_combine(dvd, weight_fac, vp0, normal, poly_r, d_r, coef_lo, coef_n, mapped, q_right)
                        #:endif

                        #:if hllc_mp_weno and WENO_RADIUS == 2
                            curvature_lo = vp0 + vm2 - vm1*2._wp
                            curvature_mid = vp1 + vm1 - vp0*2._wp
                            curvature_hi = vp2 + vp0 - vp1*2._wp
                            #:for PAIR, A, B in (('lo', 'curvature_lo', 'curvature_mid'), ('hi', 'curvature_mid', 'curvature_hi'))
                                four_a_b = 4._wp*${A}$ - ${B}$
                                four_b_a = 4._wp*${B}$ - ${A}$
                                curvature_sign = sign(1._wp, four_a_b)
                                limited_${PAIR}$ = (curvature_sign + sign(1._wp, four_b_a))*abs((curvature_sign + sign(1._wp, &
                                                    & ${A}$))*(curvature_sign + sign(1._wp, ${B}$)))*min(abs(four_a_b), &
                                                    & abs(${A}$), abs(four_b_a), abs(${B}$))/8._wp
                            #:endfor
                            #:for SIDE, SIGN, DELTA, NEIGHBOR, MD, LC in (('left', '-', 'vp1 - vp0', 'vm1', 'lo', 'hi'), &
                                                  & ('right', '+', 'vp0 - vm1', 'vp1', 'hi', 'lo'))
                                bound_ul = vp0 ${SIGN}$ (${DELTA}$)*2._wp
                                bound_md = (vp0 + ${NEIGHBOR}$ - limited_${MD}$)*5.e-1_wp
                                bound_lc = vp0 ${SIGN}$ (${DELTA}$)*5.e-1_wp + (4._wp/3._wp)*limited_${LC}$
                                bound_min = max(min(vp0, ${NEIGHBOR}$, bound_md), min(vp0, bound_ul, bound_lc))
                                bound_max = min(max(vp0, ${NEIGHBOR}$, bound_md), max(vp0, bound_ul, bound_lc))
                                q_${SIDE}$ = min(max(q_${SIDE}$, bound_min), bound_max)
                            #:endfor
                        #:endif
                    #:endif

                end subroutine s_fused_weno_state
            #:endif
