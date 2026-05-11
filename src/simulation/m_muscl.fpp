!>
!! @file
!! @brief Contains module m_muscl

#:include 'macros.fpp'

!> @brief MUSCL reconstruction with interface sharpening for contact-preserving advection
module m_muscl

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
#ifdef MFC_OpenACC
    use openacc
#endif

    use m_mpi_proxy
    use m_helper
    use m_thinc
    use m_nvtx

    private; public :: s_initialize_muscl_module, s_muscl, s_finalize_muscl_module

    integer :: v_size
    $:GPU_DECLARE(create='[v_size]')

    type(int_bounds_info) :: is1_muscl, is2_muscl, is3_muscl
    $:GPU_DECLARE(create='[is1_muscl, is2_muscl, is3_muscl]')

    !> @name The cell-average variables that will be MUSCL-reconstructed. Formerly, they are stored in v_vf. However, they are
    !! transferred to v_rs_wsL and v_rs_wsR as to be reshaped (RS) and/or characteristically decomposed. The reshaping allows the
    !! muscl procedure to be independent of the coordinate direction of the reconstruction. Lastly, notice that the left (L) and
    !! right (R) results of the characteristic decomposition are stored in custom-constructed muscl- stencils (WS) that are annexed
    !! to each position of a given scalar field.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: v_rs_ws_x_muscl, v_rs_ws_y_muscl, v_rs_ws_z_muscl
    !> @}
    $:GPU_DECLARE(create='[v_rs_ws_x_muscl, v_rs_ws_y_muscl, v_rs_ws_z_muscl]')

contains

    !> Allocate and initialize MUSCL reconstruction working arrays
    subroutine s_initialize_muscl_module()

        ! Initializing in x-direction
        is1_muscl%beg = -buff_size; is1_muscl%end = m - is1_muscl%beg
        if (n == 0) then
            is2_muscl%beg = 0
        else
            is2_muscl%beg = -buff_size
        end if

        is2_muscl%end = n - is2_muscl%beg

        if (p == 0) then
            is3_muscl%beg = 0
        else
            is3_muscl%beg = -buff_size
        end if

        is3_muscl%end = p - is3_muscl%beg

        @:ALLOCATE(v_rs_ws_x_muscl(is1_muscl%beg:is1_muscl%end, is2_muscl%beg:is2_muscl%end, is3_muscl%beg:is3_muscl%end, &
                   & 1:sys_size))

        if (n == 0) return

        ! initializing in y-direction
        is2_muscl%beg = -buff_size; is2_muscl%end = n - is2_muscl%beg
        is1_muscl%beg = -buff_size; is1_muscl%end = m - is1_muscl%beg

        if (p == 0) then
            is3_muscl%beg = 0
        else
            is3_muscl%beg = -buff_size
        end if

        is3_muscl%end = p - is3_muscl%beg

        @:ALLOCATE(v_rs_ws_y_muscl(is2_muscl%beg:is2_muscl%end, is1_muscl%beg:is1_muscl%end, is3_muscl%beg:is3_muscl%end, &
                   & 1:sys_size))

        if (p == 0) return

        ! initializing in z-direction
        is2_muscl%beg = -buff_size; is2_muscl%end = n - is2_muscl%beg
        is1_muscl%beg = -buff_size; is1_muscl%end = m - is1_muscl%beg
        is3_muscl%beg = -buff_size; is3_muscl%end = p - is3_muscl%beg

        @:ALLOCATE(v_rs_ws_z_muscl(is3_muscl%beg:is3_muscl%end, is2_muscl%beg:is2_muscl%end, is1_muscl%beg:is1_muscl%end, &
                   & 1:sys_size))

    end subroutine s_initialize_muscl_module

    !> Perform MUSCL reconstruction of left and right cell-boundary values from cell-averaged variables
    subroutine s_muscl(v_vf, vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, muscl_dir, is1_muscl_d, &

        & is2_muscl_d, is3_muscl_d)

        type(scalar_field), dimension(1:), intent(in)                                          :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_rs_vf_x, vL_rs_vf_y, &
             & vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z
        integer, intent(in)               :: muscl_dir
        type(int_bounds_info), intent(in) :: is1_muscl_d, is2_muscl_d, is3_muscl_d
        integer                           :: j, k, l, i
        real(wp)                          :: slopeL, slopeR, slope

        is1_muscl = is1_muscl_d
        is2_muscl = is2_muscl_d
        is3_muscl = is3_muscl_d

        $:GPU_UPDATE(device='[is1_muscl, is2_muscl, is3_muscl]')

        if (muscl_order /= 1 .or. dummy) then
            call s_initialize_muscl(v_vf, muscl_dir)
        end if

        if (muscl_order == 1 .or. dummy) then
            if (muscl_dir == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, ubound(v_vf, 1)
                    do l = is3_muscl%beg, is3_muscl%end
                        do k = is2_muscl%beg, is2_muscl%end
                            do j = is1_muscl%beg, is1_muscl%end
                                vL_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                vR_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (muscl_dir == 2) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                do i = 1, ubound(v_vf, 1)
                    do l = is3_muscl%beg, is3_muscl%end
                        do k = is2_muscl%beg, is2_muscl%end
                            do j = is1_muscl%beg, is1_muscl%end
                                vL_rs_vf_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                vR_rs_vf_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (muscl_dir == 3) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                do i = 1, ubound(v_vf, 1)
                    do l = is3_muscl%beg, is3_muscl%end
                        do k = is2_muscl%beg, is2_muscl%end
                            do j = is1_muscl%beg, is1_muscl%end
                                vL_rs_vf_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                vR_rs_vf_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        if (muscl_order == 2 .or. dummy) then
            ! MUSCL Reconstruction
            #:for MUSCL_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (muscl_dir == ${MUSCL_DIR}$) then
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[i, j, k, l, slopeL, slopeR, slope]')
                    do l = is3_muscl%beg, is3_muscl%end
                        do k = is2_muscl%beg, is2_muscl%end
                            do j = is1_muscl%beg, is1_muscl%end
                                do i = 1, v_size
                                    slopeL = v_rs_ws_${XYZ}$_muscl(j + 1, k, l, i) - v_rs_ws_${XYZ}$_muscl(j, k, l, i)
                                    slopeR = v_rs_ws_${XYZ}$_muscl(j, k, l, i) - v_rs_ws_${XYZ}$_muscl(j - 1, k, l, i)
                                    slope = 0._wp

                                    if (muscl_lim == 0) then  ! unlimited (central difference)
                                        slope = 5e-1_wp*(slopeL + slopeR)
                                    else if (muscl_lim == 1) then  ! minmod
                                        if (slopeL*slopeR > muscl_eps) then
                                            slope = min(abs(slopeL), abs(slopeR))
                                        end if
                                        if (slopeL < 0._wp) slope = -slope
                                    else if (muscl_lim == 2) then  ! MC
                                        if (slopeL*slopeR > muscl_eps) then
                                            slope = min(2._wp*abs(slopeL), 2._wp*abs(slopeR))
                                            slope = min(slope, 5e-1_wp*(abs(slopeL) + abs(slopeR)))
                                        end if
                                        if (slopeL < 0._wp) slope = -slope
                                    else if (muscl_lim == 3) then  ! Van Albada
                                        if (slopeL*slopeR > muscl_eps) then
                                            slope = ((slopeL + slopeR)*slopeL*slopeR)/(slopeL**2._wp + slopeR**2._wp)
                                        end if
                                    else if (muscl_lim == 4) then  ! Van Leer
                                        if (slopeL*slopeR > muscl_eps) then
                                            slope = 2._wp*slopeL*slopeR/(slopeL + slopeR)
                                        end if
                                    else if (muscl_lim == 5) then  ! SUPERBEE
                                        if (slopeL*slopeR > muscl_eps) then
                                            slope = -1._wp*min(-min(2._wp*abs(slopeL), abs(slopeR)), -min(abs(slopeL), &
                                                               & 2._wp*abs(slopeR)))
                                        end if
                                    end if

                                    ! reconstruct from left side
                                    vL_rs_vf_${XYZ}$ (j, k, l, i) = v_rs_ws_${XYZ}$_muscl(j, k, l, i) - (5.e-1_wp*slope)

                                    ! reconstruct from the right side
                                    vR_rs_vf_${XYZ}$ (j, k, l, i) = v_rs_ws_${XYZ}$_muscl(j, k, l, i) + (5.e-1_wp*slope)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endfor
        end if

        if (int_comp > 0 .and. v_size >= eqn_idx%adv%end) then
            call nvtxStartRange("WENO-INTCOMP")
            #:for MUSCL_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (muscl_dir == ${MUSCL_DIR}$) then
                    call s_thinc_compression(v_rs_ws_${XYZ}$_muscl, vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, &
                                             & vR_rs_vf_z, muscl_dir, is1_muscl, is2_muscl, is3_muscl)
                end if
            #:endfor
            call nvtxEndRange()
        end if

    end subroutine s_muscl

    !> Reshape cell-averaged variable data into direction-local work arrays for MUSCL reconstruction
    subroutine s_initialize_muscl(v_vf, muscl_dir)

        type(scalar_field), dimension(:), intent(in) :: v_vf
        integer, intent(in)                          :: muscl_dir
        integer                                      :: j, k, l, q  !< Generic loop iterators
        ! Determine MUSCL-reconstructed variables and map coordinate directions

        v_size = ubound(v_vf, 1)
        $:GPU_UPDATE(device='[v_size]')

        if (muscl_dir == 1) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=4)
            do j = 1, v_size
                do q = is3_muscl%beg, is3_muscl%end
                    do l = is2_muscl%beg, is2_muscl%end
                        do k = is1_muscl%beg - muscl_polyn, is1_muscl%end + muscl_polyn
                            v_rs_ws_x_muscl(k, l, q, j) = v_vf(j)%sf(k, l, q)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Reshaping/Projecting onto Characteristic Fields in y-direction
        if (n == 0) return

        if (muscl_dir == 2) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=4)
            do j = 1, v_size
                do q = is3_muscl%beg, is3_muscl%end
                    do l = is2_muscl%beg, is2_muscl%end
                        do k = is1_muscl%beg - muscl_polyn, is1_muscl%end + muscl_polyn
                            v_rs_ws_y_muscl(k, l, q, j) = v_vf(j)%sf(l, k, q)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Reshaping/Projecting onto Characteristic Fields in z-direction
        if (p == 0) return
        if (muscl_dir == 3) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=4)
            do j = 1, v_size
                do q = is3_muscl%beg, is3_muscl%end
                    do l = is2_muscl%beg, is2_muscl%end
                        do k = is1_muscl%beg - muscl_polyn, is1_muscl%end + muscl_polyn
                            v_rs_ws_z_muscl(k, l, q, j) = v_vf(j)%sf(q, l, k)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_initialize_muscl

    !> Finalize the MUSCL module
    subroutine s_finalize_muscl_module()

        @:DEALLOCATE(v_rs_ws_x_muscl)

        if (n == 0) return

        @:DEALLOCATE(v_rs_ws_y_muscl)

        if (p == 0) return

        @:DEALLOCATE(v_rs_ws_z_muscl)

    end subroutine s_finalize_muscl_module

end module m_muscl
