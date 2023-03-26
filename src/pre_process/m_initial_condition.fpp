!>
!! @file m_initial_condition.f90
!! @brief Contains module m_initial_condition

!> @brief This module provides a platform that is analagous to constructive
!!              solid geometry techniques and in this way allows for the creation
!!              of a wide variety of initial conditions. Several 1D, 2D and 3D
!!              fundamental geometries are included that may further be combined
!!              into more complex shapes. This is achieved by carefully setting
!!              up the order in which the patches are laid out in the domain and
!!              specifying the priority that each patch has over the preceeding
!!              ones. The resulting shapes may be identified both by the values
!!              of their primitive variables and the associated patch identities.
!!              Note that the user may choose to read in and modify a preexisting
!!              initial condition. The module m_start_up.f90 is responsible for
!!             reading in the relevant data files.
module m_initial_condition

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_patches

    use m_assign_variables
    ! ==========================================================================
    ! ==========================================================================

    implicit none

    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !< primitive variables

    type(scalar_field), allocatable, dimension(:) :: q_cons_vf !< conservative variables

    integer, allocatable, dimension(:, :, :) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.

contains

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    subroutine s_initialize_initial_condition_module() ! -------------------

        integer :: i !< generic loop iterator

        ! Allocating the primitive and conservative variables
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(0:m, 0:n, 0:p))
            allocate (q_cons_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n, 0:p))

        ! Setting default values for conservative and primitive variables so
        ! that in the case that the initial condition is wrongly laid out on
        ! the grid the simulation component will catch the problem on start-
        ! up. The conservative variables do not need to be similarly treated
        ! since they are computed directly from the primitive variables.
        do i = 1, sys_size
            q_cons_vf(i)%sf = dflt_real
            q_prim_vf(i)%sf = dflt_real
        end do

        ! Setting default values for patch identities bookkeeping variable.
        ! This is necessary to avoid any confusion in the assessment of the
        ! extent of application that the overwrite permissions give a patch
        ! when it is being applied in the domain.
        patch_id_fp = 0

    end subroutine s_initialize_initial_condition_module ! -----------------

    !>  This subroutine peruses the patches and depending on the
        !!              type of geometry associated with a particular patch, it
        !!              calls the related subroutine to setup the said geometry
        !!              on the grid using the primitive variables included with
        !!              the patch parameters. The subroutine is complete once the
        !!              primitive variables are converted to conservative ones.
    subroutine s_generate_initial_condition() ! ----------------------------

        integer :: i  !< Generic loop operator
		
        ! Converting the conservative variables to the primitive ones given
        ! preexisting initial condition data files were read in on start-up
        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                               q_prim_vf)
        end if

        !  3D Patch Geometries =============================================
        if (p > 0) then

            do i = 1, num_patches

                ! Spherical patch
                if (patch_icpp(i)%geometry == 8) then
                    call s_sphere(i, patch_id_fp, q_prim_vf)

                    ! Cuboidal patch
                elseif (patch_icpp(i)%geometry == 9) then
                    call s_cuboid(i, patch_id_fp, q_prim_vf)

                    ! Cylindrical patch
                elseif (patch_icpp(i)%geometry == 10) then
                    call s_cylinder(i, patch_id_fp, q_prim_vf)

                    ! Swept plane patch
                elseif (patch_icpp(i)%geometry == 11) then
                    call s_sweep_plane(i, patch_id_fp, q_prim_vf)

                    ! Ellipsoidal patch
                elseif (patch_icpp(i)%geometry == 12) then
                    call s_ellipsoid(i, patch_id_fp, q_prim_vf)

                    ! Analytical function patch for testing purposes
                elseif (patch_icpp(i)%geometry == 13) then
                    call s_3D_analytical(i, patch_id_fp, q_prim_vf)

                    ! Spherical harmonic patch
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_spherical_harmonic(i, patch_id_fp, q_prim_vf)

                    ! 3D Modified circular patch
                elseif (patch_icpp(i)%geometry == 19) then
                    call s_3dvarcircle(i, patch_id_fp, q_prim_vf)

                end if

            end do

            ! ==================================================================

            ! 2D Patch Geometries ==============================================
        elseif (n > 0) then

            do i = 1, num_patches

                ! Circular patch
                if (patch_icpp(i)%geometry == 2) then
                    call s_circle(i, patch_id_fp, q_prim_vf)

                    ! Rectangular patch
                elseif (patch_icpp(i)%geometry == 3) then
                    call s_rectangle(i, patch_id_fp, q_prim_vf)

                    ! Swept line patch
                elseif (patch_icpp(i)%geometry == 4) then
                    call s_sweep_line(i, patch_id_fp, q_prim_vf)

                    ! Elliptical patch
                elseif (patch_icpp(i)%geometry == 5) then
                    call s_ellipse(i, patch_id_fp, q_prim_vf)

                    ! Isentropic vortex patch
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_isentropic_vortex(i, patch_id_fp, q_prim_vf)

                    ! Analytical function patch for testing purposes
                elseif (patch_icpp(i)%geometry == 7) then
                    call s_2D_analytical(i, patch_id_fp, q_prim_vf)

                    ! Spiral patch
                elseif (patch_icpp(i)%geometry == 17) then
                    call s_spiral(i, patch_id_fp, q_prim_vf)

                    ! Modified circular patch
                elseif (patch_icpp(i)%geometry == 18) then
                    call s_varcircle(i, patch_id_fp, q_prim_vf)

                end if

            end do

            ! ==================================================================

            ! 1D Patch Geometries ==============================================
        else

            do i = 1, num_patches

                ! Line segment patch
                if (patch_icpp(i)%geometry == 1) then
                    call s_line_segment(i, patch_id_fp, q_prim_vf)

                    ! 1d analytical
                elseif (patch_icpp(i)%geometry == 15) then
                    call s_1d_analytical(i, patch_id_fp, q_prim_vf)

                    ! 1d bubble screen with sinusoidal pressure pulse
                elseif (patch_icpp(i)%geometry == 16) then
                    call s_1d_bubble_pulse(i, patch_id_fp, q_prim_vf)
                end if

            end do

        end if
        ! ==================================================================

        if (perturb_flow) call s_perturb_surrounding_flow()
        if (perturb_sph) call s_perturb_sphere()

		if (instability_wave) call s_superposition_instability_wave()

        ! Converting the primitive variables to the conservative ones
        call s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                           q_cons_vf)

    end subroutine s_generate_initial_condition ! --------------------------

    subroutine s_perturb_sphere() ! ----------------------------------------

        integer :: i, j, k, l !< generic loop operators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: alpha_unadv
        real(kind(0d0)) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
!                        IF ((perturb_alpha >= 25d-2) .AND. (perturb_alpha <= 75d-2)) THEN
                    if ((perturb_alpha /= 0d0) .and. (perturb_alpha /= 1d0)) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere ! --------------------------------------

    subroutine s_perturb_surrounding_flow() ! ------------------------------

        integer :: i, j, k, l !<  generic loop iterators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    perturb_alpha = q_prim_vf(E_idx + perturb_flow_fluid)%sf(i, j, k)
                    ! IF (perturb_alpha == 1d0) THEN
                    ! Perturb partial density
!                            CALL RANDOM_NUMBER(rand_real)
!                            rand_real = rand_real / 1d2 / 1d3
!                            q_prim_vf(perturb_flow_fluid)%sf(i,j,k) = q_prim_vf(perturb_flow_fluid)%sf(i,j,k) + rand_real
                    ! Perturb velocity
                    call random_number(rand_real)
                    rand_real = rand_real*1.d-2
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    q_prim_vf(mom_idx%end)%sf(i, j, k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    if (bubbles) then
                        q_prim_vf(alf_idx)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(alf_idx)%sf(i, j, k)
                    end if
                    ! END IF
                end do
            end do
        end do

    end subroutine s_perturb_surrounding_flow ! ----------------------------
	
	subroutine s_superposition_instability_wave() ! ------------------------------
		real(kind(0d0)), dimension(5,0:m,0:n,0:p) :: wave,wave1,wave2,wave_tmp
		real(kind(0d0)) :: tr,ti
		real(kind(0d0)) :: Ly
		integer :: i,j,k

		Ly = y_domain%end - y_domain%beg
		
		write(*,*) "generate instability waves ..."

		wave = 0d0
		wave1 = 0d0
		wave2 = 0d0
		if (p .eq. 0) then
			call s_instability_wave(2*pi*4.0/Ly,0d0,tr,ti,wave_tmp,0d0)
			wave = wave + wave_tmp
			call s_instability_wave(2*pi*2.0/Ly,0d0,tr,ti,wave_tmp,0d0)
			wave = wave + wave_tmp
			call s_instability_wave(2*pi*1.0/Ly,0d0,tr,ti,wave_tmp,0d0)
			wave = wave + wave_tmp
			wave = wave*0.05
		else
			call s_instability_wave(2*pi*4.0/Ly,0d0,tr,ti,wave_tmp,0d0)
			wave1 = wave1 + wave_tmp
			call s_instability_wave(2*pi*2.0/Ly,0d0,tr,ti,wave_tmp,0d0)
			wave1 = wave1 + wave_tmp
			call s_instability_wave(2*pi*1.0/Ly,0d0,tr,ti,wave_tmp,0d0)
			wave1 = wave1 + wave_tmp
			
			call s_instability_wave(2*pi*4.0/Ly, 2*pi*4.0/Ly,tr,ti,wave_tmp,2*pi*0.8147d0)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*2.0/Ly, 2*pi*2.0/Ly,tr,ti,wave_tmp,2*pi*0.9058d0)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*1.0/Ly, 2*pi*1.0/Ly,tr,ti,wave_tmp,2*pi*0.0270d0)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*4.0/Ly,-2*pi*4.0/Ly,tr,ti,wave_tmp,2*pi*0.5324d0)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*2.0/Ly,-2*pi*2.0/Ly,tr,ti,wave_tmp,2*pi*0.2975d0)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*1.0/Ly,-2*pi*1.0/Ly,tr,ti,wave_tmp,2*pi*0.7634d0)
			wave2 = wave2 + wave_tmp

			wave = 0.05*wave1+0.15*wave2
		end if

        do k = 0, p
        do j = 0, n
        do i = 0, m
			q_prim_vf(mom_idx%beg  )%sf(i,j,k) = q_prim_vf(mom_idx%beg  )%sf(i,j,k)+wave(2,i,j,k)
			q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = q_prim_vf(mom_idx%beg+1)%sf(i,j,k)+wave(3,i,j,k)
			if (p .gt. 0) then
			q_prim_vf(mom_idx%beg+2)%sf(i,j,k) = q_prim_vf(mom_idx%beg+2)%sf(i,j,k)+wave(4,i,j,k)
			end if
		end do
		end do
		end do

	end subroutine s_superposition_instability_wave

    subroutine s_instability_wave(alpha,beta,tr,ti,wave,shift)
		real(kind(0d0)),intent(in) :: alpha, beta !<  spatial wavenumbers
		real(kind(0d0)),dimension(0:n) :: rho_mean, u_mean, t_mean !<  mean profiles
		real(kind(0d0)),dimension(0:n) :: drho_mean, du_mean, dt_mean !<  mean profiles
		real(kind(0d0)),dimension(0:n,0:n) :: d
		real(kind(0d0)),dimension(0:5*(n+1)-1,0:5*(n+1)-1) :: ar,ai,br,bi,ci
		real(kind(0d0)),dimension(0:5*(n+1)-1,0:5*(n+1)-1) :: zr,zi
		real(kind(0d0)),dimension(0:5*(n+1)-1) :: wr,wi,fv1,fv2,fv3
		real(kind(0d0)) :: tr,ti
		real(kind(0d0)),dimension(0:5*(n+1)-1) :: vr,vi,vnr,vni
		real(kind(0d0)),dimension(5,0:m,0:n,0:p) :: wave
		real(kind(0d0)) :: shift
		
		integer :: ierr
		integer :: j, k, l !<  generic loop iterators
		integer :: ii, jj !< block matrix indicies

		real(kind(0d0)) :: gam,pi_inf,rho1,mach,c1

		gam = 1.+1./fluid_pp(1)%gamma
		pi_inf = fluid_pp(1)%pi_inf*(gam-1.)/gam
		rho1 = patch_icpp(1)%alpha_rho(1)/patch_icpp(1)%alpha(1)
		c1 = sqrt((gam*(patch_icpp(1)%pres+pi_inf))/rho1)
		mach = 1./c1

		! Assign mean profiles
		do j=0,n
			u_mean(j)=tanh(y_cc(j))
			t_mean(j)=1+0.5*(gam-1)*mach**2*(1-u_mean(j)**2)
			rho_mean(j)=1/T_mean(j)
		end do

		! Compute differential operator
		dy = y_cc(1)-y_cc(0)
		d=0d0
		d(1,0)=-1/(2*dy)
		d(1,2)= 1/(2*dy)
		do j=2,n-2
			d(j,j-2)= 1/(12*dy)
			d(j,j-1)=-8/(12*dy)
			d(j,j+1)= 8/(12*dy)
			d(j,j+2)=-1/(12*dy)
		end do
		d(n-1,n-2)=-1/(2*dy)
		d(n-1,n)  = 1/(2*dy)

		! Compute y-derivatives of rho, u, T
		do j=0,n
			drho_mean(j)=0
			du_mean(j)=0
			dt_mean(j)=0
			do k=0,n
				drho_mean(j) = drho_mean(j)+d(j,k)*rho_mean(k)
				du_mean(j) = du_mean(j)+d(j,k)*u_mean(k)
				dt_mean(j) = dt_mean(j)+d(j,k)*t_mean(k)
			end do
		end do

		! Compute B and C -> A=B+C
		br=0d0
		bi=0d0
		ci=0d0
		do j=0,n
		    ii = 1; jj = 1; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*u_mean(j); 
			ii = 1; jj = 2; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*rho_mean(j);
			ii = 1; jj = 3; bi((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = -drho_mean(j);
			ii = 1; jj = 4; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = beta*rho_mean(j);

			ii = 2; jj = 1; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*t_mean(j)/(rho_mean(j)*gam*mach**2);
			ii = 2; jj = 2; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*u_mean(j);
			ii = 2; jj = 3; bi((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = -du_mean(j);
			ii = 2; jj = 5; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha/(gam*mach**2);

			ii = 3; jj = 1; bi((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = -dt_mean(j)/(rho_mean(j)*gam*mach**2);
			ii = 3; jj = 3; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*u_mean(j);
			ii = 3; jj = 5; bi((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = -drho_mean(j)/(rho_mean(j)*gam*mach**2);

			ii = 4; jj = 1; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = beta*t_mean(j)/(rho_mean(j)*gam*mach**2);
			ii = 4; jj = 4; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*u_mean(j);
			ii = 4; jj = 5; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = beta/(gam*mach**2);

			ii = 5; jj = 2; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = (gam-1)*alpha/rho_mean(j);
			ii = 5; jj = 3; bi((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = -dt_mean(j);
			ii = 5; jj = 4; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = (gam-1)*beta/rho_mean(j);
			ii = 5; jj = 5; br((ii-1)*(n+1)+j,(jj-1)*(n+1)+j) = alpha*u_mean(j);

			do k=0,n
				ii = 1; jj = 3; ci((ii-1)*(n+1)+j,(jj-1)*(n+1)+k) = -rho_mean(j)*d(j,k);
				ii = 3; jj = 1; ci((ii-1)*(n+1)+j,(jj-1)*(n+1)+k) = -t_mean(j)*d(j,k)/(rho_mean(j)*gam*mach**2);
				ii = 3; jj = 5; ci((ii-1)*(n+1)+j,(jj-1)*(n+1)+k) = -d(j,k)/(gam*mach**2);
				ii = 5; jj = 3; ci((ii-1)*(n+1)+j,(jj-1)*(n+1)+k) = -(gam-1)*d(j,k)/rho_mean(j);
			end do
		end do
		ar = br
		ai = bi+ci

		! Compute Eigenvalues & Eigenfunctions
		call cg(5*(n+1),5*(n+1),ar,ai,wr,wi,zr,zi,fv1,fv2,fv3,ierr)

		! Find the most unstable wave
		call find_unstable_mode(5*(n+1),wr,wi,zr,zi,tr,ti,vr,vi)

		! Normalize
		call normalize_eigvec(5*(n+1),vr,vi,vnr,vni)

		! Generate instability waves
		call generate_wave(5*(n+1),vnr,vni,alpha,beta,wave,shift)

    end subroutine s_instability_wave

	subroutine generate_wave(nl,vnr,vni,alpha,beta,wave,shift)
		integer nl
		real(kind(0d0)) :: alpha,beta,ang
		real(kind(0d0)), dimension(0:nl-1) :: vnr,vni
		real(kind(0d0)), dimension(5,0:m,0:n,0:p) :: wave
		real(kind(0d0)) :: shift
		integer i,j,k
		
		do i=0,m
		do j=0,n
		do k=0,p
			if (beta .eq. 0) then
				ang = alpha*x_cc(i)
			else 
				ang = alpha*x_cc(i)+beta*z_cc(k)+shift
			end if
			
			wave(1,i,j,k) = vnr(j)*cos(ang) &
						   -vni(j)*sin(ang)		! rho
			wave(2,i,j,k) = vnr((n+1)+j)*cos(ang) &
						   -vni((n+1)+j)*sin(ang)	! u
			wave(3,i,j,k) = vnr(2*(n+1)+j)*cos(ang) &
						   -vni(2*(n+1)+j)*sin(ang)	! v
			wave(4,i,j,k) = vnr(3*(n+1)+j)*cos(ang) &
						   -vni(3*(n+1)+j)*sin(ang)	! w
			wave(5,i,j,k) = vnr(4*(n+1)+j)*cos(ang) &
						   -vni(4*(n+1)+j)*sin(ang)	! T
		end do
		end do
		end do
	
	end subroutine generate_wave
	
	subroutine normalize_eigvec(nl,vr,vi,vnr,vni)
	integer nl
	real(kind(0d0)), dimension(nl) :: vr,vi,vnr,vni
	real(kind(0d0)) :: norm
	real(kind(0d0)) :: tr,ti,cr,ci
	integer i,idx
	
	norm = 0d0
	do i=1,nl
		if (dsqrt(vr(i)**2+vi(i)**2) .gt. dabs(norm)) then
			idx = i
			if (vr(i) .gt. 0) then 
				norm = sqrt(vr(i)**2+vi(i)**2)
			else
				norm =-sqrt(vr(i)**2+vi(i)**2)
			end if
		end if
	end do

	vnr = vr/norm
	vni = vi/norm

	tr = vnr(idx)
	ti = vni(idx)
	do i=1,nl
		call cdiv(vnr(i),vni(i),tr,ti,cr,ci)
		vnr(i) = cr
		vni(i) = ci
	end do
	
	end subroutine normalize_eigvec

	subroutine find_unstable_mode(nl,wr,wi,zr,zi,tr,ti,vr,vi)
		integer nl
		real(kind(0d0)), dimension(nl) :: wr,wi
		real(kind(0d0)), dimension(nl,nl) :: zr,zi
		real(kind(0d0)) :: tr,ti
		real(kind(0d0)), dimension(nl) :: vr,vi
		integer i,k
		
		k=1
		do i=2,nl
			if (wi(i) .gt. wi(k)) then
				k = i
			end if
		end do
	
		tr = wr(k)
		ti = wi(k)
		vr = zr(:,k)
		vi = zi(:,k)
	
	end subroutine find_unstable_mode

!====================================================================
!	Subroutines & functions for computing eigenvalues and eigenvectors
!	which are modified from EISPACK (https://netlib.org/eispack/)
!		1) cg
!		2) cbal
!		3) corth
!		4) comqr2
!		5) csroot
!		6) cdiv
!		7) pythag
!====================================================================
	subroutine cg(nm,nl,ar,ai,wr,wi,zr,zi,fv1,fv2,fv3,ierr)
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        nl  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex general matrix.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for comqr
!           and comqr2.  the normal completion code is zero.
!
!        fv1, fv2, and  fv3  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------	
		integer nm,nl,is1,is2,ierr
		real(kind(0d0)), dimension(nm,nl) :: ar,ai,zr,zi
		real(kind(0d0)), dimension(nl) :: wr,wi,fv1,fv2,fv3

		if (nl .le. nm) go to 10
		ierr = 10*nl
		go to 50

10		call cbal(nm,nl,ar,ai,is1,is2,fv1)
		call corth(nm,nl,is1,is2,ar,ai,fv2,fv3)
		call comqr2(nm,nl,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
		if (ierr .ne. 0) go to 50
		call cbabk2(nm,nl,is1,is2,fv1,nl,zr,zi)

50	return
    end subroutine cg

    subroutine cbal(nm,nl,ar,ai,low,igh,scale)
!     this subroutine is a translation of the algol procedure
!     cbalance, which is a complex version of balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a complex matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        nl is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex matrix to be balanced.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the balanced matrix.
!
!        low and igh are two integers such that ar(i,j) and ai(i,j)
!          are equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,nl.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j)       j = low,...,igh
!                 = p(j)         j = igh+1,...,nl.
!     the order in which the interchanges are made is nl to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in cbalance appears in
!     cbal  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     arithmetic is real throughout.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
		integer i,j,k,l,ml,nl,jj,nm,igh,low,iexc
		real(kind(0d0)), dimension(nm,nl) :: ar,ai
		real(kind(0d0)), dimension(nl) :: scale
		real(kind(0d0)) :: c,f,g,r,s,b2,radix
		logical noconv

		radix = 16.0d0

		b2 = radix * radix
		k = 1
		l = nl
		go to 100
!     .......... in-line procedure for row and
!                column exchange ..........
20  	scale(ml) = j
		if (j .eq. ml) go to 50

		do 30 i = 1, l
			f = ar(i,j)
			ar(i,j) = ar(i,ml)
			ar(i,ml) = f
			f = ai(i,j)
			ai(i,j) = ai(i,ml)
			ai(i,ml) = f
30  	continue

		do 40 i = k, nl
			f = ar(j,i)
			ar(j,i) = ar(ml,i)
			ar(ml,i) = f
			f = ai(j,i)
			ai(j,i) = ai(ml,i)
			ai(ml,i) = f
40  	continue

50  	go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
80  	if (l .eq. 1) go to 280
			l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
100 	do 120 jj = 1, l
			j = l + 1 - jj

        do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
110     continue

        ml = l
        iexc = 1
        go to 20
120 	continue

		go to 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
130 	k = k + 1

140 	do 170 j = k, l

        do 150 i = k, l
           if (i .eq. j) go to 150
           if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
150     continue

        ml = k
        iexc = 2
        go to 20
170 	continue
!     .......... now balance the submatrix in rows k to l ..........
		do 180 i = k, l
		scale(i) = 1.0d0
180		continue
!     .......... iterative loop for norm reduction ..........
190 	noconv = .false.

		do 270 i = k, l
			c = 0.0d0
			r = 0.0d0

        do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
200     continue
!     .......... guard against zero c or r due to underflow ..........
        if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
        g = r / radix
        f = 1.0d0
        s = c + r
210     if (c .ge. g) go to 220
        f = f * radix
        c = c * b2
        go to 210
220     g = r * radix
230     if (c .lt. g) go to 240
        f = f / radix
        c = c / b2
        go to 230
!     .......... now balance ..........
240     if ((c + r) / f .ge. 0.95d0 * s) go to 270
        g = 1.0d0 / f
        scale(i) = scale(i) * f
        noconv = .true.

        do 250 j = k, nl
           ar(i,j) = ar(i,j) * g
           ai(i,j) = ai(i,j) * g
250     continue

        do 260 j = 1, l
           ar(j,i) = ar(j,i) * f
           ai(j,i) = ai(j,i) * f
260     continue

270 	continue

		if (noconv) go to 190

280 	low = k
		igh = l
	return
    end subroutine cbal
	
    subroutine corth(nm,nl,low,igh,ar,ai,ortr,orti)
!     this subroutine is a translation of a complex analogue of
!     the algol procedure orthes, num. math. 12, 349-368(1968)
!     by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a complex general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     unitary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        nl is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=nl.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex input matrix.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the hessenberg matrix.  information
!          about the unitary transformations used in the reduction
!          is stored in the remaining triangles under the
!          hessenberg matrix.
!
!        ortr and orti contain further information about the
!          transformations.  only elements low through igh are used.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
		integer i,j,ml,nl,ii,jj,la,mp,nm,igh,kp1,low
		real(kind(0d0)),dimension(nm,nl) :: ar,ai
		real(kind(0d0)),dimension(igh) :: ortr,orti
		real(kind(0d0)) :: f,g,h,fi,fr,scale,c

		integer mll
		mll = 6
		
		la = igh - 1
		kp1 = low + 1
		if (la .lt. kp1) go to 200

		do 180 ml = kp1, la
			h = 0.0d0
			ortr(ml) = 0.0d0
			orti(ml) = 0.0d0
			scale = 0.0d0
!     .......... scale column (algol tol then not needed) ..........
        do 90 i = ml, igh
        scale = scale + dabs(ar(i,ml-1)) + dabs(ai(i,ml-1))
90		continue
        if (scale .eq. 0d0) go to 180
        mp = ml + igh
!     .......... for i=igh step -1 until ml do -- ..........
        do 100 ii = ml, igh
            i = mp - ii
            ortr(i) = ar(i,ml-1) / scale
            orti(i) = ai(i,ml-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
100     continue
!
        g = dsqrt(h)
		call pythag(ortr(ml),orti(ml),f)
        if (f .eq. 0d0) go to 103
        h = h + f * g
        g = g / f
        ortr(ml) = (1.0d0 + g) * ortr(ml)
        orti(ml) = (1.0d0 + g) * orti(ml)
        go to 105

103     ortr(ml) = g
        ar(ml,ml-1) = scale
!     .......... form (i-(u*ut)/h) * a ..........
105     do 130 j = ml, nl
            fr = 0.0d0
            fi = 0.0d0
!     .......... for i=igh step -1 until ml do -- ..........
            do 110 ii = ml, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
110         continue
!
            fr = fr / h
            fi = fi / h
!
            do 120 i = ml, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
120         continue
!
130     continue
!     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
        do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
!     .......... for j=igh step -1 until ml do -- ..........
            do 140 jj = ml, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
140         continue
!
            fr = fr / h
            fi = fi / h
!
            do 150 j = ml, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
150         continue
!
160     continue
!
        ortr(ml) = scale * ortr(ml)
        orti(ml) = scale * orti(ml)
        ar(ml,ml-1) = -g * ar(ml,ml-1)
        ai(ml,ml-1) = -g * ai(ml,ml-1)
180     continue
!
200 return
    end subroutine corth
	
    subroutine comqr2(nm,nl,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
!  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
!  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
!
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a complex upper hessenberg matrix by the qr
!     method.  the eigenvectors of a complex general matrix
!     can also be found if  corth  has been used to reduce
!     this general matrix to hessenberg form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        nl is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=nl.
!
!        ortr and orti contain information about the unitary trans-
!          formations used in the reduction by  corth, if performed.
!          only elements low through igh are used.  if the eigenvectors
!          of the hessenberg matrix are desired, set ortr(j) and
!          orti(j) to 0.0d0 for these elements.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain further
!          information about the transformations which were used in the
!          reduction by  corth, if performed.  if the eigenvectors of
!          the hessenberg matrix are desired, these elements may be
!          arbitrary.
!
!     on output
!
!        ortr, orti, and the upper hessenberg portions of hr and hi
!          have been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,nl.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors.  the eigenvectors
!          are unnormalized.  if an error exit is made, none of
!          the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*nl iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated october 1989.
!
!     ------------------------------------------------------------------
		integer i,j,k,l,ml,nl,en,ii,jj,ll,nm,nn,igh,ip1,&
              itn,its,low,lp1,enm1,iend,ierr
		real(kind(0d0)),dimension(nm,nl) :: hr,hi,zr,zi
		real(kind(0d0)),dimension(nl) :: wr,wi
		real(kind(0d0)),dimension(igh) :: ortr,orti
		real(kind(0d0)) :: si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,&
						norm,tst1,tst2,c,d
!
		ierr = 0
!     .......... initialize eigenvector matrix ..........
		do 101 j = 1, nl
!
        do 100 i = 1, nl
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
100   	continue
        zr(j,j) = 1.0d0
101 	continue
!     .......... form the matrix of accumulated transformations
!                from the information left by corth ..........
		iend = igh - low - 1
		if (iend .lt. 0) go to 180
		if (iend .eq. 0) go to 150
		if (iend .gt. 0) go to 105
!     .......... for i=igh-1 step -1 until low+1 do -- ..........
105 	do 140 ii = 1, iend
        i = igh - ii
        if (dabs(ortr(i)) .eq. 0d0 .and. dabs(orti(i)) .eq. 0d0) go to 140
        if (dabs(hr(i,i-1)) .eq. 0d0 .and. dabs(hi(i,i-1)) .eq. 0d0) go to 140
!     .......... norm below is negative of h formed in corth ..........
        norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
        ip1 = i + 1

        do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
110   	continue
!
        do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
!
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
115         continue
!
            sr = sr / norm
            si = si / norm
!
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
120         continue
!
130     continue
!
140     continue
!     .......... create real subdiagonal elements ..........
150     l = low + 1
!
        do 170 i = l, igh
        ll = min0(i+1,igh)
        if (dabs(hi(i,i-1)) .eq. 0d0) go to 170
        call pythag(hr(i,i-1),hi(i,i-1),norm)
        yr = hr(i,i-1) / norm
        yi = hi(i,i-1) / norm
        hr(i,i-1) = norm
        hi(i,i-1) = 0.0d0
!
        do 155 j = i, nl
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
155     continue
!
        do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
160     continue
!
        do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
165     continue

170     continue
!     .......... store roots isolated by cbal ..........
180     do 200 i = 1, nl
        if (i .ge. low .and. i .le. igh) go to 200
        wr(i) = hr(i,i)
        wi(i) = hi(i,i)
200     continue
!
		en = igh
		tr = 0.0d0
		ti = 0.0d0
		itn = 30*nl
!     .......... search for next eigenvalue ..........
220 	if (en .lt. low) go to 680
		its = 0
		enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
240 	do 260 ll = low, en
        l = en + low - ll
        if (l .eq. low) go to 300
        tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1)) &
             + dabs(hr(l,l)) + dabs(hi(l,l))
        tst2 = tst1 + dabs(hr(l,l-1))
        if (tst2 .eq. tst1) go to 300
260 	continue
!     .......... form shift ..........
300 	if (l .eq. en) go to 660
		if (itn .eq. 0) go to 1000
		if (its .eq. 10 .or. its .eq. 20) go to 320
		sr = hr(en,en)
		si = hi(en,en)
		xr = hr(enm1,en) * hr(en,enm1)
		xi = hi(enm1,en) * hr(en,enm1)
		if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
		yr = (hr(enm1,enm1) - sr) / 2.0d0
		yi = (hi(enm1,enm1) - si) / 2.0d0
		call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
		if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
		zzr = -zzr
		zzi = -zzi
310 	call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
		sr = sr - xr
		si = si - xi
		go to 340
!     .......... form exceptional shift ..........
320 	sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
		si = 0.0d0
!
340 	do 360 i = low, en
        hr(i,i) = hr(i,i) - sr
        hi(i,i) = hi(i,i) - si
360 	continue
!
		tr = tr + sr
		ti = ti + si
		its = its + 1
		itn = itn - 1
!     .......... reduce to triangle (rows) ..........
		lp1 = l + 1
!
		do 500 i = lp1, en
        sr = hr(i,i-1)
        hr(i,i-1) = 0.0d0
		call pythag(hr(i-1,i-1),hi(i-1,i-1),c)
		call pythag(c,sr,norm)
        xr = hr(i-1,i-1) / norm
        wr(i-1) = xr
        xi = hi(i-1,i-1) / norm
        wi(i-1) = xi
        hr(i-1,i-1) = norm
        hi(i-1,i-1) = 0.0d0
        hi(i,i-1) = sr / norm
!
        do 490 j = i, nl
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
490     continue
!
500 	continue
!
		si = hi(en,en)
		if (dabs(si) .eq. 0d0) go to 540
		call pythag(hr(en,en),si,norm)
		sr = hr(en,en) / norm
		si = si / norm
		hr(en,en) = norm
		hi(en,en) = 0.0d0
		if (en .eq. nl) go to 540
		ip1 = en + 1
!
		do 520 j = ip1, nl
        yr = hr(en,j)
        yi = hi(en,j)
        hr(en,j) = sr * yr + si * yi
        hi(en,j) = sr * yi - si * yr
520 	continue
!     .......... inverse operation (columns) ..........
540 	do 600 j = lp1, en
        xr = wr(j-1)
        xi = wi(j-1)
!
        do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
560         hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
580     continue
!
        do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
590     continue

600 	continue
!
        if (dabs(si) .eq. 0d0) go to 240
!
        do 630 i = 1, en
        yr = hr(i,en)
        yi = hi(i,en)
        hr(i,en) = sr * yr - si * yi
        hi(i,en) = sr * yi + si * yr
630 	continue
!
		do 640 i = low, igh
        yr = zr(i,en)
        yi = zi(i,en)
        zr(i,en) = sr * yr - si * yi
        zi(i,en) = sr * yi + si * yr
640 	continue
!
		go to 240
!     .......... a root found ..........
660 	hr(en,en) = hr(en,en) + tr
		wr(en) = hr(en,en)
		hi(en,en) = hi(en,en) + ti
		wi(en) = hi(en,en)
		en = enm1
		go to 220
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
680 	norm = 0.0d0
!
		do i = 1, nl
        do j = i, nl
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
		end do
		end do
!
		if (nl .eq. 1 .or. norm .eq. 0d0) go to 1001
!     .......... for en=nl step -1 until 2 do -- ..........
		do 800 nn = 2, nl
        en = nl + 2 - nn
        xr = wr(en)
        xi = wi(en)
        hr(en,en) = 1.0d0
        hi(en,en) = 0.0d0
        enm1 = en - 1
!     .......... for i=en-1 step -1 until 1 do -- ..........
        do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1

            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
740         continue
!
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
760            yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
765         continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
!     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
770         continue
!
780	    continue
!
800		continue
!     .......... end backsubstitution ..........
!     .......... vectors of isolated roots ..........
		do  840 i = 1, nl
        if (i .ge. low .and. i .le. igh) go to 840
!
        do 820 j = I, nl
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
820     continue
!
840 	continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=nl step -1 until low do -- ..........
		do jj = low, nl
        j = nl + low - jj
        ml = min0(j,igh)
!
        do i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
!
            do 860 k = low, ml
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
860         continue
!
            zr(i,j) = zzr
            zi(i,j) = zzi
		end do
		end do
!
		go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*nl iterations ..........
1000 	ierr = en
1001 	return
    end subroutine comqr2

	subroutine cbabk2(nm,nl,low,igh,scale,ml,zr,zi)
!
      integer i,j,k,ml,nl,ii,nm,igh,low
      double precision scale(nl),zr(nm,ml),zi(nm,ml)
      double precision s
!
!     this subroutine is a translation of the algol procedure
!     cbabk2, which is a complex version of balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  cbal.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        nl is the order of the matrix.
!
!        low and igh are integers determined by  cbal.
!
!        scale contains information determining the permutations
!          and scaling factors used by  cbal.
!
!        ml is the number of eigenvectors to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first ml columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first ml columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (ml .eq. 0) go to 200
      if (igh .eq. low) go to 120
!
      do 110 i = low, igh
         s = scale(i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0d0/scale(i). ..........
         do 100 j = 1, ml
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
100    continue
!
110 continue
!     .......... for i=low-1 step -1 until 1,
!                igh+1 step 1 until nl do -- ..........
120 do 140 ii = 1, nl
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
!
         do 130 j = 1, ml
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
130    continue
!
140 continue
!
200 return
	end subroutine cbabk2
	  
    subroutine csroot(xr,xi,yr,yi)
		real(kind(0d0)) :: xr,xi,yr,yi
!
!     (yr,yi) = complex dsqrt(xr,xi) 
!     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
!
		real(kind(0d0)) :: s,tr,ti,c
		tr = xr
		ti = xi
		call pythag(tr,ti,c)
		s = dsqrt(0.5d0*(c + dabs(tr)))
		if (tr .ge. 0.0d0) yr = s
		if (ti .lt. 0.0d0) s = -s
		if (tr .le. 0.0d0) yi = s
		if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
		if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
		return
    end subroutine csroot
	  
	subroutine cdiv(ar,ai,br,bi,cr,ci)
		real(kind(0d0)) :: ar,ai,br,bi,cr,ci
!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
		real(kind(0d0)) :: s,ars,ais,brs,bis
		s = dabs(br) + dabs(bi)
		ars = ar/s
		ais = ai/s
		brs = br/s
		bis = bi/s
		s = brs**2 + bis**2
		cr = (ars*brs + ais*bis)/s
		ci = (ais*brs - ars*bis)/s
		return
    end subroutine cdiv

	subroutine pythag(a,b,c)
		real(kind(0d0)) :: a,b,c
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
		real(kind(0d0)) :: p,r,s,t,u
		p = dmax1(dabs(a),dabs(b))
		if (p .eq. 0.0d0) go to 20
		r = (dmin1(dabs(a),dabs(b))/p)**2
10 		continue
        t = 4.0d0 + r
        if (t .eq. 4.0d0) go to 20
        s = r/t
        u = 1.0d0 + 2.0d0*s
        p = u*p
        r = (s/u)**2 * r
        go to 10
20 		c = p
		return
    end subroutine pythag




    !>  Deallocation procedures for the module
    subroutine s_finalize_initial_condition_module() ! ---------------------

        integer :: i !< Generic loop iterator

        ! Dellocating the primitive and conservative variables
        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        ! Deallocating the patch identities bookkeeping variable
        deallocate (patch_id_fp)

    end subroutine s_finalize_initial_condition_module ! -------------------

end module m_initial_condition
