!>
!! @file m_fftw.f90
!! @brief Contains module m_fftw

!> @brief The module contains the subroutines for the FFT routines
MODULE m_fftw

    ! Dependencies =============================================================
    USE, INTRINSIC :: ISO_C_BINDING

    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy

#IFDEF _OPENACC
    USE cudafor

    USE cufft
#ENDIF

    ! ==========================================================================

    IMPLICIT NONE

    PRIVATE; PUBLIC :: s_initialize_fftw_module,   &
                       s_apply_fourier_filter,     &
                       s_finalize_fftw_module
    
    INCLUDE 'fftw3.f03'

    TYPE(C_PTR) :: fwd_plan, bwd_plan
    TYPE(C_PTR) :: fftw_real_data, fftw_cmplx_data, fftw_fltr_cmplx_data
    INTEGER :: real_size, cmplx_size, x_size, batch_size

    REAL(C_DOUBLE), POINTER :: data_real(:) !< Real data

    COMPLEX(C_DOUBLE_COMPLEX), POINTER :: data_cmplx(:) !< 
    !! Complex data in Fourier space
    
    

    COMPLEX(C_DOUBLE_COMPLEX), POINTER :: data_fltr_cmplx(:) !< 
    !! Filtered complex data in Fourier space



    
    
#IFDEF _OPENACC
    REAL(kind(0d0)), POINTER :: data_real_gpu(:) 

    COMPLEX(kind(0d0)), POINTER :: data_cmplx_gpu(:)  

    COMPLEX(kind(0d0)), POINTER :: data_fltr_cmplx_gpu(:) 

    INTEGER :: fwd_plan_gpu, bwd_plan_gpu, ierr
     
    INTEGER, allocatable :: cufft_size(:), iembed(:), oembed(:)

    INTEGER :: istride, ostride, idist, odist, rank
     
    !$acc declare create(real_size, cmplx_size, x_size, batch_size, data_real_gpu, data_cmplx_gpu, data_fltr_cmplx_gpu)
#ENDIF

    CONTAINS




        !>  The purpose of this subroutine is to create the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
        SUBROUTINE s_initialize_fftw_module() ! ----------------------------------



            ! Size of input array going into DFT
            real_size = p+1
            ! Size of output array coming out of DFT
            cmplx_size = (p+1)/2+1

            x_size = m+1

            batch_size = x_size*sys_size

#IFDEF _OPENACC
            rank = 1; istride = 1; ostride = 1

            allocate(cufft_size(1:rank),iembed(1:rank), oembed(1:rank))

            cufft_size(1) = real_size; 

            iembed(1) = 0
            oembed(1) = 0

            !$acc update device(real_size, cmplx_size, x_size, sys_size, batch_size)
#ENDIF
            ! Allocate input and output DFT data sizes
            fftw_real_data       = fftw_alloc_real   (int( real_size, C_SIZE_T))
            fftw_cmplx_data      = fftw_alloc_complex(int(cmplx_size, C_SIZE_T))
            fftw_fltr_cmplx_data = fftw_alloc_complex(int(cmplx_size, C_SIZE_T))
            ! Associate input and output data pointers with allocated memory
            CALL c_f_pointer(fftw_real_data ,      data_real ,      [ real_size])
            CALL c_f_pointer(fftw_cmplx_data,      data_cmplx,      [cmplx_size])
            CALL c_f_pointer(fftw_fltr_cmplx_data, data_fltr_cmplx, [cmplx_size])

            ! Generate plans for forward and backward DFTs
            fwd_plan = fftw_plan_dft_r2c_1d(real_size, data_real      , data_cmplx, FFTW_ESTIMATE)
            bwd_plan = fftw_plan_dft_c2r_1d(real_size, data_fltr_cmplx, data_real , FFTW_ESTIMATE)

#IFDEF _OPENACC
            allocate(data_real_gpu(1:real_size*x_size*sys_size))
            allocate(data_cmplx_gpu(1:cmplx_size*x_size*sys_size))
            allocate(data_fltr_cmplx_gpu(1:cmplx_size*x_size*sys_size))

            ierr = cufftPlanMany(fwd_plan_gpu, rank, cufft_size, iembed, istride, real_size, oembed, ostride, cmplx_size, CUFFT_D2Z, batch_size)
            ierr = cufftPlanMany(bwd_plan_gpu, rank, cufft_size, iembed, istride, cmplx_size, oembed, ostride, real_size, CUFFT_Z2D, batch_size)
#ENDIF    
        END SUBROUTINE s_initialize_fftw_module ! ------------------------------




        !>  The purpose of this subroutine is to apply a Fourier low-
        !!      pass filter to the flow variables in the azimuthal direction
        !!      to remove the high-frequency content. This alleviates the 
        !!      restrictive CFL condition arising from cells near the axis.
        SUBROUTINE s_apply_fourier_filter(q_cons_vf) ! --------------------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf

            INTEGER :: Nfq !< Number of kept modes

            INTEGER :: i,j,k,l !< Generic loop iterators

            ! Restrict filter to processors that have cells adjacent to axis
            IF (bc_y%beg >= 0) RETURN

#IFDEF _OPENACC
            
!$acc parallel loop collapse(3) gang vector default(present) 
                DO k = 1, sys_size
                    DO j = 0, m
                        DO l = 1, cmplx_size
                            data_fltr_cmplx_gpu(l + j*cmplx_size + (k-1)*cmplx_size*x_size) = (0d0, 0d0)
                        END DO
                    END DO
                END DO


!$acc parallel loop collapse(3) gang vector default(present) 
                DO k = 1, sys_size
                    DO j = 0, m
                        DO l = 0, p
                            data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size) = q_cons_vf(k)%sf(j, 0, l)
                        END DO
                    END DO
                END DO

!$acc host_data use_device(data_real_gpu, data_cmplx_gpu)                
                ierr = cufftExecD2Z(fwd_plan_gpu, data_real_gpu, data_cmplx_gpu)
!$acc end host_data

                Nfq = 3

!$acc parallel loop collapse(3) gang vector default(present) firstprivate(Nfq) 
                DO k = 1, sys_size
                     DO j = 0, m
                        DO l = 1, Nfq
                            data_fltr_cmplx_gpu(l + j*cmplx_size + (k-1)*cmplx_size*x_size) = data_cmplx_gpu(l + j*cmplx_size+ (k-1)*cmplx_size*x_size)
                        END DO 
                    END DO
                END DO

!$acc host_data use_device(data_real_gpu, data_fltr_cmplx_gpu)  
                ierr = cufftExecZ2D(bwd_plan_gpu, data_fltr_cmplx_gpu, data_real_gpu)
!$acc end host_data

                

!$acc parallel loop collapse(3) gang vector default(present)
                DO k = 1, sys_size
                    DO j = 0, m 
                        DO l = 0, p
                            data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size) = data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size)/REAL(real_size,KIND(0d0))
                            q_cons_vf(k)%sf(j, 0, l) = data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size)
                        END DO
                    END DO
                END DO
           
            
            DO i = 1, fourier_rings

!$acc parallel loop collapse(3) gang vector default(present) 
                DO k = 1, sys_size
                    DO j = 0, m
                        DO l = 1, cmplx_size
                            data_fltr_cmplx_gpu(l + j*cmplx_size + (k-1)*cmplx_size*x_size) = (0d0, 0d0)
                        END DO
                    END DO
                END DO

!$acc parallel loop collapse(3) gang vector default(present) firstprivate(i)
                DO k = 1, sys_size
                    DO j = 0, m
                        DO l = 0, p
                            data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size) = q_cons_vf(k)%sf(j, i, l)
                        END DO
                    END DO
                END DO

!$acc host_data use_device(data_real_gpu, data_cmplx_gpu)  
                ierr = cufftExecD2Z(fwd_plan_gpu, data_real_gpu, data_cmplx_gpu)
!$acc end host_data

                Nfq = MIN(FLOOR(2d0*REAL(i,KIND(0d0))*pi),cmplx_size)

!$acc parallel loop collapse(3) gang vector default(present) firstprivate(Nfq)
                DO k = 1, sys_size 
                    DO j = 0, m
                        DO l = 1, Nfq
                            data_fltr_cmplx_gpu(l + j*cmplx_size + (k-1)*cmplx_size*x_size) = data_cmplx_gpu(l + j*cmplx_size + (k-1)*cmplx_size*x_size)
                        END DO
                    END DO 
                END DO

!$acc host_data use_device(data_real_gpu, data_fltr_cmplx_gpu)  
                ierr = cufftExecZ2D(bwd_plan_gpu, data_fltr_cmplx_gpu, data_real_gpu)
!$acc end host_data

!$acc parallel loop collapse(3) gang vector default(present) firstprivate(i)
                DO k = 1, sys_size
                    DO j = 0, m
                        DO l = 0, p
                            data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size) = data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size)/REAL(real_size,KIND(0d0))
                            q_cons_vf(k)%sf(j, i, l) = data_real_gpu(l + j*real_size + 1 + (k-1)*real_size*x_size)
                        END DO
                    END DO 
                END DO

    
            END DO
            
              
#ELSE
            Nfq = 3
            DO j = 0, m
                DO k = 1, sys_size
                    data_fltr_cmplx(:) = (0d0,0d0)
                    data_real(1:p+1) = q_cons_vf(k)%sf(j,0,0:p)
                    CALL fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                    data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                    CALL fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                    data_real(:) = data_real(:)/REAL(real_size,KIND(0d0))
                    q_cons_vf(k)%sf(j,0,0:p) = data_real(1:p+1)
                END DO
            END DO

            ! Apply Fourier filter to additional rings
            DO i = 1, fourier_rings
                Nfq = MIN(FLOOR(2d0*REAL(i,KIND(0d0))*pi),cmplx_size)
                DO j = 0, m
                    DO k = 1, sys_size
                        data_fltr_cmplx(:) = (0d0,0d0)
                        data_real(1:p+1) = q_cons_vf(k)%sf(j,i,0:p)
                        CALL fftw_execute_dft_r2c(fwd_plan, data_real, data_cmplx)
                        data_fltr_cmplx(1:Nfq) = data_cmplx(1:Nfq)
                        CALL fftw_execute_dft_c2r(bwd_plan, data_fltr_cmplx, data_real)
                        data_real(:) = data_real(:)/REAL(real_size,KIND(0d0))
                        q_cons_vf(k)%sf(j,i,0:p) = data_real(1:p+1)
                    END DO
                END DO
            END DO              
#ENDIF



        END SUBROUTINE s_apply_fourier_filter ! --------------------------------




        !>  The purpose of this subroutine is to destroy the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
        SUBROUTINE s_finalize_fftw_module() ! ------------------------------------
            CALL fftw_free(fftw_real_data)
            CALL fftw_free(fftw_cmplx_data)
            CALL fftw_free(fftw_fltr_cmplx_data)

            CALL fftw_destroy_plan(fwd_plan)
            CALL fftw_destroy_plan(bwd_plan)

#IFDEF _OPENACC
            deallocate(data_real_gpu, data_fltr_cmplx_gpu, data_cmplx_gpu)
            ierr = cufftDestroy(fwd_plan_gpu)
            ierr = cufftDestroy(bwd_plan_gpu)
#ENDIF    

        END SUBROUTINE s_finalize_fftw_module ! --------------------------------



END MODULE
