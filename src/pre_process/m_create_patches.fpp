module m_create_patches

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_assign_patches
    ! ==========================================================================

    implicit none


    real(kind(0d0)) :: radius
    real(kind(0d0)) :: smooth_coeff !<
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    real(kind(0d0)) :: length_x, length_y, length_z

    real(kind(0d0)) :: a, b, c, d !<
    !! When a line or a plane sweep patch geometry is employed, these variables
    !! represent the coefficients associated with the equation describing the
    !! said line or plane.

    real(kind(0d0)) :: cart_y, cart_z
    real(kind(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates
    
contains

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)

        real(kind(0d0)), intent(IN) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord ! --------------

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)

        real(kind(0d0)), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord ! --------------

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

!====================================================================
	subroutine s_superposition_instability_wave() ! ------------------------------
		real(kind(0d0)), dimension(5,0:m,0:n,0:p) :: wave,wave1,wave2,wave_tmp
		real(kind(0d0)) :: tr,ti
		integer :: i,j,k

		write(*,*) "generate instability waves ..."

		wave = 0d0
		wave1 = 0d0
		wave2 = 0d0
		if (p .eq. 0) then
			call s_instability_wave(2*pi*4.0/59.0,0d0,tr,ti,wave_tmp,0d0)
			wave = wave + wave_tmp
			call s_instability_wave(2*pi*2.0/59.0,0d0,tr,ti,wave_tmp,0d0)
			wave = wave + wave_tmp
			call s_instability_wave(2*pi*1.0/59.0,0d0,tr,ti,wave_tmp,0d0)
			wave = wave + wave_tmp
			wave = wave*0.05
		else
			call s_instability_wave(2*pi*4.0/59.0,0d0,tr,ti,wave_tmp,0d0)
			wave1 = wave1 + wave_tmp
			call s_instability_wave(2*pi*2.0/59.0,0d0,tr,ti,wave_tmp,0d0)
			wave1 = wave1 + wave_tmp
			call s_instability_wave(2*pi*1.0/59.0,0d0,tr,ti,wave_tmp,0d0)
			wave1 = wave1 + wave_tmp
			
			call s_instability_wave(2*pi*4.0/59.0, 2*pi*4.0/59.0,tr,ti,wave_tmp,11d0/31d0*2*pi)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*2.0/59.0, 2*pi*2.0/59.0,tr,ti,wave_tmp,13d0/31d0*2*pi)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*1.0/59.0, 2*pi*1.0/59.0,tr,ti,wave_tmp,17d0/31d0*2*pi)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*4.0/59.0,-2*pi*4.0/59.0,tr,ti,wave_tmp,19d0/31d0*2*pi)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*2.0/59.0,-2*pi*2.0/59.0,tr,ti,wave_tmp,23d0/31d0*2*pi)
			wave2 = wave2 + wave_tmp
			call s_instability_wave(2*pi*1.0/59.0,-2*pi*1.0/59.0,tr,ti,wave_tmp,29d0/31d0*2*pi)
			wave2 = wave2 + wave_tmp
			! call s_instability_wave(2*pi*4.0/59.0, 2*pi*4.0/59.0,tr,ti,wave_tmp,0.8147d0*2*pi)
			! wave2 = wave2 + wave_tmp
			! call s_instability_wave(2*pi*2.0/59.0, 2*pi*2.0/59.0,tr,ti,wave_tmp,0.9058d0*2*pi)
			! wave2 = wave2 + wave_tmp
			! call s_instability_wave(2*pi*1.0/59.0, 2*pi*1.0/59.0,tr,ti,wave_tmp,0.0270d0*2*pi)
			! wave2 = wave2 + wave_tmp
			! call s_instability_wave(2*pi*4.0/59.0,-2*pi*4.0/59.0,tr,ti,wave_tmp,0.5324d0*2*pi)
			! wave2 = wave2 + wave_tmp
			! call s_instability_wave(2*pi*2.0/59.0,-2*pi*2.0/59.0,tr,ti,wave_tmp,0.2975d0*2*pi)
			! wave2 = wave2 + wave_tmp
			! call s_instability_wave(2*pi*1.0/59.0,-2*pi*1.0/59.0,tr,ti,wave_tmp,0.7634d0*2*pi)
			! wave2 = wave2 + wave_tmp

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

!====================================================================
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

		! Find the most unstable wave and normalize it
		call find_unstable_mode(5*(n+1),wr,wi,zr,zi,tr,ti,vr,vi)

		! Normalize
		call normalize_eigvec(5*(n+1),vr,vi,vnr,vni)

		! Generate instability waves
		call generate_wave(5*(n+1),vnr,vni,alpha,beta,wave,shift)

    end subroutine s_instability_wave ! ----------------------------

!====================================================================
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
	
!====================================================================
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

!====================================================================
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

!====================================================================
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
	
!====================================================================
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
	
!====================================================================
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

!====================================================================
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
	  
!====================================================================
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
	  
!====================================================================
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

!====================================================================
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

    !>          The line segment patch is a 1D geometry that may be used,
        !!              for example, in creating a Riemann problem. The geometry
        !!              of the patch is well-defined when its centroid and length
        !!              in the x-coordinate direction are provided. Note that the
        !!              line segment patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !! @param patch_id patch identifier
    subroutine s_line_segment(patch_id) ! ----------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma

        integer :: i, j  !< Generic loop operators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the line segment's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and end x-coordinates of the line segment
        ! based on its centroid and length
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the line segment patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0)

                !IF ( (q_prim_vf(1)%sf(i,0,0) < 1.e-12) .AND. (model_eqns .NE. 4)) THEN
                !    !zero density, reassign according to Tait EOS
                !    q_prim_vf(1)%sf(i,0,0) = &
                !        (((q_prim_vf(E_idx)%sf(i,0,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                !        rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,0,0))
                !END IF
            end if
        end do

    end subroutine s_line_segment ! ----------------------------------------

    !>  The spiral patch is a 2D geometry that may be used, The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !!  @param patch_id patch identifier
    subroutine s_spiral(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< Generic loop iterators
        real(kind(0d0)) :: th, thickness, nturns, mya
        real(kind(0d0)) :: spiral_x_min, spiral_x_max, spiral_y_min, spiral_y_max

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        mya = patch_icpp(patch_id)%radius
        thickness = patch_icpp(patch_id)%length_x
        nturns = patch_icpp(patch_id)%length_y

!
        logic_grid = 0
        do k = 0, int(m*91*nturns)
            th = k/real(int(m*91d0*nturns))*nturns*2.d0*pi

            spiral_x_min = minval((/f_r(th, 0.0d0, mya)*cos(th), &
                                    f_r(th, thickness, mya)*cos(th)/))
            spiral_y_min = minval((/f_r(th, 0.0d0, mya)*sin(th), &
                                    f_r(th, thickness, mya)*sin(th)/))

            spiral_x_max = maxval((/f_r(th, 0.0d0, mya)*cos(th), &
                                    f_r(th, thickness, mya)*cos(th)/))
            spiral_y_max = maxval((/f_r(th, 0.0d0, mya)*sin(th), &
                                    f_r(th, thickness, mya)*sin(th)/))

            do j = 0, n; do i = 0, m; 
                    if ((x_cc(i) > spiral_x_min) .and. (x_cc(i) < spiral_x_max) .and. &
                        (y_cc(j) > spiral_y_min) .and. (y_cc(j) < spiral_y_max)) then
                        logic_grid(i, j, 0) = 1
                    end if
                end do; end do
        end do

        do j = 0, n
            do i = 0, m
                if ((logic_grid(i, j, 0) == 1)) then
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)
                end if
            end do
        end do

    end subroutine s_spiral ! ----------------------------------------------

    !> Archimedes spiral function
        !! @param myth Angle
        !! @param offset Thickness
        !! @param a Starting position
    function f_r(myth, offset, a)
        real(kind(0d0)), intent(IN) :: myth, offset, a
        real(kind(0d0)) :: b
        real(kind(0d0)) :: f_r

        !r(th) = a + b*th

        b = 2.d0*a/(2.d0*pi)
        f_r = a + b*myth + offset
    end function f_r

    !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_circle(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j !< Generic loop iterators

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then

                    eta = tanh(smooth_coeff/min(dx, dy)* &
                               (sqrt((x_cc(i) - x_centroid)**2 &
                                     + (y_cc(j) - y_centroid)**2) &
                                - radius))*(-0.5d0) + 0.5d0

                end if

                if (((x_cc(i) - x_centroid)**2 &
                     + (y_cc(j) - y_centroid)**2 <= radius**2 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                end if

            end do
        end do

    end subroutine s_circle ! ----------------------------------------------

    !>             The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
    subroutine s_varcircle(patch_id) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j

        real(kind(0d0)) :: myr, thickness

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        thickness = patch_icpp(patch_id)%epsilon

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                myr = dsqrt((x_cc(i) - x_centroid)**2 &
                            + (y_cc(j) - y_centroid)**2)

                if (myr <= radius + thickness/2.d0 .and. &
                    myr >= radius - thickness/2.d0 .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                    q_prim_vf(alf_idx)%sf(i, j, 0) = patch_icpp(patch_id)%alpha(1)* &
                                                     dexp(-0.5d0*((myr - radius)**2.d0)/(thickness/3.d0)**2.d0)
                end if

            end do
        end do

    end subroutine s_varcircle ! ----------------------------------------------

    subroutine s_3dvarcircle(patch_id) ! ----------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j, k

        real(kind(0d0)) :: myr, thickness

        ! Transferring the circular patch's radius, centroid, smearing patch
        ! identity and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_z = patch_icpp(patch_id)%length_z
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff
        thickness = patch_icpp(patch_id)%epsilon

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the circular patch's boundary is enabled.
        eta = 1d0

        ! write for all z

        ! Checking whether the circle covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    myr = dsqrt((x_cc(i) - x_centroid)**2 &
                                + (y_cc(j) - y_centroid)**2)

                    if (myr <= radius + thickness/2.d0 .and. &
                        myr >= radius - thickness/2.d0 .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                        q_prim_vf(alf_idx)%sf(i, j, k) = patch_icpp(patch_id)%alpha(1)* &
                                                         dexp(-0.5d0*((myr - radius)**2.d0)/(thickness/3.d0)**2.d0)
                    end if

                end do
            end do
        end do

    end subroutine s_3dvarcircle ! ----------------------------------------------

    !>      The elliptical patch is a 2D geometry. The geometry of
        !!      the patch is well-defined when its centroid and radii
        !!      are provided. Note that the elliptical patch DOES allow
        !!      for the smoothing of its boundary
        !! @param patch_id is the patch identifier
    subroutine s_ellipse(patch_id) ! ---------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j !< Generic loop operators

        ! Transferring the elliptical patch's radii, centroid, smearing
        ! patch identity, and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        a = patch_icpp(patch_id)%radii(1)
        b = patch_icpp(patch_id)%radii(2)
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value
        ! be modified as the patch is laid out on the grid, but only in
        ! the case that smoothing of the elliptical patch's boundary is
        ! enabled.
        eta = 1d0

        ! Checking whether the ellipse covers a particular cell in the
        ! domain and verifying whether the current patch has permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = tanh(smooth_coeff/min(dx, dy)* &
                               (sqrt(((x_cc(i) - x_centroid)/a)**2 + &
                                     ((y_cc(j) - y_centroid)/b)**2) &
                                - 1d0))*(-0.5d0) + 0.5d0
                end if

                if ((((x_cc(i) - x_centroid)/a)**2 + &
                     ((y_cc(j) - y_centroid)/b)**2 <= 1d0 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)
                end if
            end do
        end do

    end subroutine s_ellipse ! ---------------------------------------------

    !>      The ellipsoidal patch is a 3D geometry. The geometry of
        !!       the patch is well-defined when its centroid and radii
        !!       are provided. Note that the ellipsoidal patch DOES allow
        !!       for the smoothing of its boundary
        !! @param patch_id is the patch identifier
    subroutine s_ellipsoid(patch_id) ! -------------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j, k

        ! Transferring the ellipsoidal patch's radii, centroid, smearing
        ! patch identity, and smearing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        a = patch_icpp(patch_id)%radii(1)
        b = patch_icpp(patch_id)%radii(2)
        c = patch_icpp(patch_id)%radii(3)
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value
        ! be modified as the patch is laid out on the grid, but only in
        ! the case that smoothing of the ellipsoidal patch's boundary is
        ! enabled.
        eta = 1d0

        ! Checking whether the ellipsoid covers a particular cell in the
        ! domain and verifying whether the current patch has permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then
                        eta = tanh(smooth_coeff/min(dx, dy, dz)* &
                                   (sqrt(((x_cc(i) - x_centroid)/a)**2 + &
                                         ((cart_y - y_centroid)/b)**2 + &
                                         ((cart_z - z_centroid)/c)**2) &
                                    - 1d0))*(-0.5d0) + 0.5d0
                    end if

                    if ((((x_cc(i) - x_centroid)/a)**2 + &
                         ((cart_y - y_centroid)/b)**2 + &
                         ((cart_z - z_centroid)/c)**2 <= 1d0 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)
                    end if
                end do
            end do
        end do

    end subroutine s_ellipsoid ! -------------------------------------------

    !>      The rectangular patch is a 2D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, in alignment with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x- and y-
        !!              coordinate directions are provided. Please note that the
        !!              rectangular patch DOES NOT allow for the smoothing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
    subroutine s_rectangle(patch_id) ! -------------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< Equation of state parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the rectangle's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates of the
        ! rectangle based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y

        ! Since the rectangular patch does not allow for its boundaries to
        ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the rectangle covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m
                if (x_boundary%beg <= x_cc(i) .and. &
                    x_boundary%end >= x_cc(i) .and. &
                    y_boundary%beg <= y_cc(j) .and. &
                    y_boundary%end >= y_cc(j) &
                    .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                    if ((q_prim_vf(1)%sf(i, j, 0) < 1.e-10) .and. (model_eqns == 4)) then
                        !zero density, reassign according to Tait EOS
                        q_prim_vf(1)%sf(i, j, 0) = &
                            (((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                            rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, j, 0))
                    end if
                end if
            end do
        end do

    end subroutine s_rectangle ! -------------------------------------------

    !>  The swept line patch is a 2D geometry that may be used,
        !!      for example, in creating a solid boundary, or pre-/post-
        !!      shock region, at an angle with respect to the axes of the
        !!      Cartesian coordinate system. The geometry of the patch is
        !!      well-defined when its centroid and normal vector, aimed
        !!      in the sweep direction, are provided. Note that the sweep
        !!      line patch DOES allow the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_sweep_line(patch_id) ! ------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j !< Generic loop operators

        ! Transferring the centroid information of the line to be swept
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Obtaining coefficients of the equation describing the sweep line
        a = patch_icpp(patch_id)%normal(1)
        b = patch_icpp(patch_id)%normal(2)
        c = -a*x_centroid - b*y_centroid

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the sweep line patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the region swept by the line covers a particular
        ! cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the
        ! primitive variables of the current patch are written to this cell.
        do j = 0, n
            do i = 0, m

                if (patch_icpp(patch_id)%smoothen) then
                    eta = 5d-1 + 5d-1*tanh(smooth_coeff/min(dx, dy) &
                                           *(a*x_cc(i) + b*y_cc(j) + c) &
                                           /sqrt(a**2 + b**2))
                end if

                if ((a*x_cc(i) + b*y_cc(j) + c >= 0d0 &
                     .and. &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    .or. &
                    patch_id_fp(i, j, 0) == smooth_patch_id) &
                    then
                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                end if

            end do
        end do

    end subroutine s_sweep_line ! ------------------------------------------

    !> The isentropic vortex is a 2D geometry that may be used,
        !!              for example, to generate an isentropic flow disturbance.
        !!              Geometry of the patch is well-defined when its centroid
        !!              and radius are provided. Notice that the patch DOES NOT
        !!              allow for the smoothing of its boundary.
        !! @param patch_id is the patch identifier
    subroutine s_isentropic_vortex(patch_id) ! ----------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j

        ! Transferring isentropic vortex patch's centroid and radius info
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        radius = patch_icpp(patch_id)%radius

        ! Since the isentropic vortex patch does not allow for its boundary
        ! to get smoothed, the pseudo volume fraction is set to 1 to ensure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Verifying whether the isentropic vortex includes a particular cell
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries work out the primitive variables of the
        ! the current patch are assigned to this cell.
        do j = 0, n
            do i = 0, m

                if ((x_cc(i) - x_centroid)**2 &
                    + (y_cc(j) - y_centroid)**2 <= radius**2 &
                    .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) &
                    then

                    call s_assign_patch_primitive_variables(patch_id, &
                                                            i, j, 0)

                end if

            end do
        end do

    end subroutine s_isentropic_vortex ! -----------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_1D_analytical(patch_id) ! ---------------------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Placeholders for the cell boundary values
        real(kind(0d0)) :: a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0)

                !what variables to alter
                !bump in pressure
                q_prim_vf(E_idx)%sf(i, 0, 0) = q_prim_vf(E_idx)%sf(i, 0, 0)* &
                                               (1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0)/(2.d0*0.005d0)))

                !bump in void fraction
                !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !bump in R(x)
                !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !IF (model_eqns == 4) THEN
                !reassign density
                !IF (num_fluids == 1) THEN
                ! q_prim_vf(1)%sf(i, 0, 0) = &
                !     (((q_prim_vf(E_idx)%sf(i, 0, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                !     rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, 0, 0))
                !END IF
                !ELSE IF (model_eqns == 2) THEN
                !can manually adjust density here
                !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )
                !END IF
            end if
        end do

    end subroutine s_1D_analytical ! ---------------------------------------

    subroutine s_1d_bubble_pulse(patch_id) ! ---------------------------------
        ! Description: This patch assigns the primitive variables as analytical
        !       functions such that the code can be verified.

        ! Patch identifier
        integer, intent(IN) :: patch_id

        ! Placeholders for the cell boundary values
        real(kind(0d0)) :: fac, a, b, c, d, pi_inf, gamma, lit_gamma

        ! Generic loop iterators
        integer :: i, j

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        length_x = patch_icpp(patch_id)%length_x

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the line segment covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do i = 0, m
            if (x_boundary%beg <= x_cc(i) .and. &
                x_boundary%end >= x_cc(i) .and. &
                patch_icpp(patch_id)%alter_patch(patch_id_fp(i, 0, 0))) then

                call s_assign_patch_primitive_variables(patch_id, i, 0, 0)

                !what variables to alter
                !sinusoid in pressure
                q_prim_vf(E_idx)%sf(i, 0, 0) = q_prim_vf(E_idx)%sf(i, 0, 0)* &
                                               (1d0 + 0.1d0*sin(-1d0*(x_cb(i) - x_centroid)*2d0*pi/length_x))

                !bump in void fraction
                !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !bump in R(x)
                !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                !IF (model_eqns == 4) THEN
                !reassign density
                !IF (num_fluids == 1) THEN
                q_prim_vf(1)%sf(i, 0, 0) = &
                    (((q_prim_vf(E_idx)%sf(i, 0, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                    rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, 0, 0))
                !END IF
                !ELSE IF (model_eqns == 2) THEN
                !can manually adjust density here
                !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )
                !END IF
            end if
        end do

    end subroutine s_1D_bubble_pulse ! ---------------------------------------

    !>  This patch assigns the primitive variables as analytical
        !!  functions such that the code can be verified.
        !!  @param patch_id is the patch identifier
    subroutine s_2D_analytical(patch_id) ! ---------------------------------

        integer, intent(IN) :: patch_id

        real(kind(0d0)) :: a, b, c, d !< placeholderrs for the cell boundary values
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters

        integer :: i, j !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y

        ! Computing the beginning and the end x- and y-coordinates
        ! of the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y

        ! Since the patch doesn't allow for its boundaries to be
        ! smoothed out, the pseudo volume fraction is set to 1 to
        ! ensure that only the current patch contributes to the fluid
        ! state in the cells that this patch covers.
        eta = 1d0

        ! Checking whether the patch covers a particular cell in the
        ! domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out,
        ! the primitive variables of the current patch are assigned
        ! to this cell.
        do j = 0, n
            do i = 0, m
                if (x_boundary%beg <= x_cc(i) .and. &
                    x_boundary%end >= x_cc(i) .and. &
                    y_boundary%beg <= y_cc(j) .and. &
                    y_boundary%end >= y_cc(j) .and. &
                    patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, 0))) then

                    call s_assign_patch_primitive_variables(patch_id, i, j, 0)

                    !what variables to alter
                    !x-y bump in pressure
                    q_prim_vf(E_idx)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)* &
                               (1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0)/(2.d0*0.005d0)))

                    !x-bump
                    !q_prim_vf(E_idx)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)* &
                    !(1d0 + 0.2d0*dexp(-1d0*((x_cb(i) - x_centroid)**2.d0)/(2.d0*0.005d0)))

                    !bump in void fraction
                    !q_prim_vf(adv_idx%beg)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !bump in R(x)
                    !q_prim_vf(adv_idx%end+1)%sf(i,j,0) = q_prim_vf(adv_idx%end+1)%sf(i,j,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !reassign density
                    !q_prim_vf(1)%sf(i, j, 0) = &
                    !(((q_prim_vf(E_idx)%sf(i, j, 0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                    !rhoref*(1d0 - q_prim_vf(alf_idx)%sf(i, j, 0))

                    ! ================================================================================

                    ! Sinusoidal initial condition for all flow variables =============================

                    ! Cell-center values
!                        a = 0d0
!                        b = 0d0
!                        c = 0d0
!                        d = 0d0
!                        q_prim_vf(adv_idx%beg)%sf(i,j,0) = SIN(x_cc(i)) * SIN(y_cc(j))
!                        q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1d0
!                        q_prim_vf(cont_idx%end)%sf(i,j,0) = (1d0 - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1d0
!                        q_prim_vf(mom_idx%beg)%sf(i,j,0) = SIN(x_cc(i))
!                        q_prim_vf(mom_idx%end)%sf(i,j,0) = SIN(y_cc(j))
!                        q_prim_vf(E_idx)%sf(i,j,0) = 1d0

                    ! Cell-average values
!                       a = x_cc(i) - 5d-1*dx ! x-beg
!                       b = x_cc(i) + 5d-1*dx ! x-end
!                       c = y_cc(j) - 5d-1*dy ! y-beg
!                       d = y_cc(j) + 5d-1*dy ! y-end
!                       q_prim_vf(adv_idx%beg)%sf(i,j,0) = 1d0/((b-a)*(d-c)) * &
!                               (COS(a)*COS(c) - COS(a)*COS(d) - COS(b)*COS(c) + COS(b)*COS(d))
!                       q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1d0
!                       q_prim_vf(cont_idx%end)%sf(i,j,0) = (1d0 - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1d0
!                       q_prim_vf(mom_idx%beg)%sf(i,j,0) = (COS(a) - COS(b))/(b-a)
!                       q_prim_vf(mom_idx%end)%sf(i,j,0) = (COS(c) - COS(d))/(d-c)
!                       q_prim_vf(E_idx)%sf(i,j,0) = 1d0
                    ! ================================================================================

                    ! Initial pressure profile smearing for bubble collapse case of Tiwari (2013) ====
                    !IF((       (x_cc(i))**2                     &
                    !         + (y_cc(j))**2 <= 1d0**2)) THEN
                    !         q_prim_vf(E_idx)%sf(i,j,0) = 1d5 / 25d0
                    !ELSE
                    !    q_prim_vf(E_idx)%sf(i,j,0) = 1d5 + 1d0/SQRT(x_cc(i)**2+y_cc(j)**2) &
                    !                                    * ((1d5/25d0) - 1d5)
                    !END IF
                    ! ================================================================================

                end if
            end do
        end do

    end subroutine s_2D_analytical ! ---------------------------------------

    !>      This patch assigns the primitive variables as analytical
        !!      functions such that the code can be verified.
        !!      @param patch_id is the patch identifier
    subroutine s_3D_analytical(patch_id) ! ---------------------------------

        integer, intent(IN) :: patch_id
        real(kind(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters

        integer :: i, j, k !< generic loop iterators

        pi_inf = fluid_pp(1)%pi_inf
        gamma = fluid_pp(1)%gamma
        lit_gamma = (1d0 + gamma)/gamma

        ! Transferring the patch's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the patch based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y
        z_boundary%beg = z_centroid - 0.5d0*length_z
        z_boundary%end = z_centroid + 0.5d0*length_z

        ! Since the analytical patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the patch covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! to that cell. If both queries check out, the primitive variables
        ! of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (x_boundary%beg <= x_cc(i) .and. &
                        x_boundary%end >= x_cc(i) .and. &
                        y_boundary%beg <= cart_y .and. &
                        y_boundary%end >= cart_y .and. &
                        z_boundary%beg <= cart_z .and. &
                        z_boundary%end >= cart_z &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                        !gaussian ball
                        !what variables to alter
                        !bump in pressure
                        q_prim_vf(E_idx)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k)* &
                                                       (1d0 + 0.2d0*exp(-1d0* &
                                                                        ((x_cb(i) - x_centroid)**2.d0 + &
                                                                         (y_cb(j) - y_centroid)**2.d0 + &
                                                                         (z_cb(k) - z_centroid)**2.d0) &
                                                                        /(2.d0*0.5d0)))

                        !bump in void fraction
                        !                       q_prim_vf(adv_idx%beg)%sf(i, j, k) = q_prim_vf(adv_idx%beg)%sf(i, j, k)* &
                        !                                                           (1d0 + 0.2d0*exp(-1d0* &
                        !                                                                           ((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0 + (z_cb(k) - z_centroid)**2.d0) &
                        !                                                                          /(2.d0*0.005d0)))

                        !bump in R(x)
                        !       q_prim_vf(adv_idx%end + 1)%sf(i, j, k) = q_prim_vf(adv_idx%end + 1)%sf(i, j, k)* &
                        !                                               (1d0 + 0.2d0*exp(-1d0* &
                        !                                                                               ((x_cb(i) - x_centroid)**2.d0 + (y_cb(j) - y_centroid)**2.d0 + (z_cb(k) - z_centroid)**2.d0) &
!                                                                                  /(2.d0*0.005d0)))

                        !reassign density
                        !         q_prim_vf(1)%sf(i, j, k) = &
                        !             (((q_prim_vf(E_idx)%sf(i, j, k) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma))* &
                        !            rhoref*(1d0 - q_prim_vf(E_idx + 1)%sf(i, j, k))

                        ! ================================================================================

                        ! Constant x-velocity in cylindrical grid ========================================
!                        q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1d0
!                        q_prim_vf(cont_idx%end )%sf(i,j,k) = 0d0
!                        q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0d0
!                        q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = COS(z_cc(k))
!                        q_prim_vf(mom_idx%end  )%sf(i,j,k) = -SIN(z_cc(k))
!                        q_prim_vf(E_idx        )%sf(i,j,k) = 1d0
!                        q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1d0
                        ! ================================================================================

                        ! Couette flow in cylindrical grid ===============================================
                        !q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1d0
                        !q_prim_vf(cont_idx%end )%sf(i,j,k) = 0d0
                        !q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0d0
                        !q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = y_cc(j)*COS(z_cc(k))*SIN(z_cc(k))
                        !q_prim_vf(mom_idx%end  )%sf(i,j,k) = -y_cc(j)*SIN(z_cc(k))**2
                        !q_prim_vf(E_idx        )%sf(i,j,k) = 1d0
                        !q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1d0
                        ! ================================================================================

                    end if

                end do
            end do
        end do

    end subroutine s_3D_analytical ! ---------------------------------------

    !>      This patch generates the shape of the spherical harmonics
        !!      as a perturbation to a perfect sphere
        !!      @param patch_id is the patch identifier
    subroutine s_spherical_harmonic(patch_id) ! ----------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< generic loop iterators

        complex(kind(0d0)) :: cmplx_i = (0d0, 1d0)
        complex(kind(0d0)) :: H

        ! Transferring the patch's centroid and radius information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        radius = patch_icpp(patch_id)%radius
        epsilon = patch_icpp(patch_id)%epsilon
        beta = patch_icpp(patch_id)%beta

        ! Since the analytical patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the patch covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! to that cell. If both queries check out, the primitive variables
        ! of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (((x_cc(i) - x_centroid)**2 &
                         + (cart_y - y_centroid)**2 &
                         + (cart_z - z_centroid)**2 <= radius**2 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k)))) &
                        then

                        call s_convert_cylindrical_to_spherical_coord(x_cc(i), y_cc(j))

                        if (epsilon == 1d0) then
                            if (beta == 0d0) then
                                H = 5d-1*sqrt(3d0/pi)*cos(sph_phi)
                            elseif (beta == 1d0) then
                                H = -5d-1*sqrt(3d0/(2d0*pi))*exp(cmplx_i*z_cc(k))*sin(sph_phi)
                            end if
                        elseif (epsilon == 2d0) then
                            if (beta == 0d0) then
                                H = 25d-2*sqrt(5d0/pi)*(3d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 1d0) then
                                H = -5d-1*sqrt(15d0/(2d0*pi))*exp(cmplx_i*z_cc(k))*sin(sph_phi)*cos(sph_phi)
                            elseif (beta == 2d0) then
                                H = 25d-2*sqrt(15d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))*sin(sph_phi)**2
                            end if
                        elseif (epsilon == 3d0) then
                            if (beta == 0d0) then
                                H = 25d-2*sqrt(7d0/pi)*(5d0*cos(sph_phi)**3d0 - 3d0*cos(sph_phi))
                            elseif (beta == 1d0) then
                                H = -125d-3*sqrt(21d0/pi)*exp(cmplx_i*z_cc(k))*sin(sph_phi)* &
                                    (5d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 2d0) then
                                H = 25d-2*sqrt(105d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*cos(sph_phi)
                            elseif (beta == 3d0) then
                                H = -125d-3*sqrt(35d0/pi)*exp(3d0*cmplx_i*z_cc(k))*sin(sph_phi)**3d0
                            end if
                        elseif (epsilon == 4d0) then
                            if (beta == 0d0) then
                                H = 3d0/16d0*sqrt(1d0/pi)*(35d0*cos(sph_phi)**4d0 - &
                                                           3d1*cos(sph_phi)**2 + 3d0)
                            elseif (beta == 1d0) then
                                H = -3d0/8d0*sqrt(5d0/pi)*exp(cmplx_i*z_cc(k))* &
                                    sin(sph_phi)*(7d0*cos(sph_phi)**3d0 - 3d0*cos(sph_phi))
                            elseif (beta == 2d0) then
                                H = 3d0/8d0*sqrt(5d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*(7d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 3d0) then
                                H = -3d0/8d0*sqrt(35d0/pi)*exp(3d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**3d0*cos(sph_phi)
                            elseif (beta == 4d0) then
                                H = 3d0/16d0*sqrt(35d0/(2d0*pi))*exp(4d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**4d0
                            end if
                        elseif (epsilon == 5d0) then
                            if (beta == 0d0) then
                                H = 1d0/16d0*sqrt(11d0/pi)*(63d0*cos(sph_phi)**5d0 - &
                                                            7d1*cos(sph_phi)**3d0 + 15d0*cos(sph_phi))
                            elseif (beta == 1d0) then
                                H = -1d0/16d0*sqrt(165d0/(2d0*pi))*exp(cmplx_i*z_cc(k))* &
                                    sin(sph_phi)*(21d0*cos(sph_phi)**4d0 - 14d0*cos(sph_phi)**2 + 1d0)
                            elseif (beta == 2d0) then
                                H = 125d-3*sqrt(1155d0/(2d0*pi))*exp(2d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**2*(3d0*cos(sph_phi)**3d0 - cos(sph_phi))
                            elseif (beta == 3d0) then
                                H = -1d0/32d0*sqrt(385d0/pi)*exp(3d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**3d0*(9d0*cos(sph_phi)**2 - 1d0)
                            elseif (beta == 4d0) then
                                H = 3d0/16d0*sqrt(385d0/(2d0*pi))*exp(4d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**4d0*cos(sph_phi)
                            elseif (beta == 5d0) then
                                H = -3d0/32d0*sqrt(77d0/pi)*exp(5d0*cmplx_i*z_cc(k))* &
                                    sin(sph_phi)**5d0
                            end if
                        end if

                        q_prim_vf(adv_idx%beg)%sf(i, j, k) = 1d0 - abs(real(H, kind(0d0)))

                    end if

                end do
            end do
        end do

    end subroutine s_spherical_harmonic ! ----------------------------------

    !>          The spherical patch is a 3D geometry that may be used,
        !!              for example, in creating a bubble or a droplet. The patch
        !!              geometry is well-defined when its centroid and radius are
        !!              provided. Please note that the spherical patch DOES allow
        !!              for the smoothing of its boundary.
        !!      @param patch_id is the patch identifier
    subroutine s_sphere(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        ! Generic loop iterators
        integer :: i, j, k !< generic loop iterators

        real(kind(0d0)) :: radius_pressure, pressure_bubble, pressure_inf !<
            !! Variables to initialize the pressure field that corresponds to the
            !! bubble-collapse test case found in Tiwari et al. (2013)

        ! Transferring spherical patch's radius, centroid, smoothing patch
        ! identity and smoothing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smoothing of the spherical patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the sphere covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! that cell. If both queries check out, the primitive variables of
        ! the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then

                        eta = tanh(smooth_coeff/min(dx, dy, dz)* &
                                   (sqrt((x_cc(i) - x_centroid)**2 &
                                         + (cart_y - y_centroid)**2 &
                                         + (cart_z - z_centroid)**2) &
                                    - radius))*(-0.5d0) + 0.5d0

                    end if

                    if (((x_cc(i) - x_centroid)**2 &
                         + (cart_y - y_centroid)**2 &
                         + (cart_z - z_centroid)**2 <= radius**2 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                    end if

                    ! Initialization of the pressure field that corresponds to the bubble-collapse
                     !! test case found in Tiwari et al. (2013)
                    ! radius_pressure = SQRT(x_cc(i)**2) ! 1D
                    ! radius_pressure = SQRT(x_cc(i)**2 + cart_y**2) ! 2D
                    ! radius_pressure = SQRT(x_cc(i)**2 + cart_y**2 + cart_z**2) ! 3D
                    ! pressure_bubble = 1.E+04
                    ! pressure_inf    = 1.E+05
                    ! q_prim_vf(E_idx)%sf(i,j,k) = pressure_inf + radius / radius_pressure * (pressure_bubble - pressure_inf)
                    !
                    ! IF(       ((  x_cc(i) - x_centroid)**2                    &
                    !          + (   cart_y - y_centroid)**2                    &
                    !          + (   cart_z - z_centroid)**2) <= radius**2)   &
                    !                               THEN
                    !
                    !    q_prim_vf(E_idx)%sf(i,j,k) = pressure_bubble
                    !
                    ! END IF

                end do
            end do
        end do

    end subroutine s_sphere ! ----------------------------------------------

    !>      The cuboidal patch is a 3D geometry that may be used, for
        !!              example, in creating a solid boundary, or pre-/post-shock
        !!              region, which is aligned with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x-, y- and
        !!              z-coordinate directions are provided. Please notice that
        !!              the cuboidal patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !!      @param patch_id is the patch identifier
    subroutine s_cuboid(patch_id) ! ----------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< Generic loop iterators

        ! Transferring the cuboid's centroid and length information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the cuboid based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y
        z_boundary%beg = z_centroid - 0.5d0*length_z
        z_boundary%end = z_centroid + 0.5d0*length_z

        ! Since the cuboidal patch does not allow for its boundaries to get
        ! smoothed out, the pseudo volume fraction is set to 1 to make sure
        ! that only the current patch contributes to the fluid state in the
        ! cells that this patch covers.
        eta = 1d0

        ! Checking whether the cuboid covers a particular cell in the domain
        ! and verifying whether the current patch has permission to write to
        ! to that cell. If both queries check out, the primitive variables
        ! of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (x_boundary%beg <= x_cc(i) .and. &
                        x_boundary%end >= x_cc(i) .and. &
                        y_boundary%beg <= cart_y .and. &
                        y_boundary%end >= cart_y .and. &
                        z_boundary%beg <= cart_z .and. &
                        z_boundary%end >= cart_z &
                        .and. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                    end if
                end do
            end do
        end do

    end subroutine s_cuboid ! ----------------------------------------------

    !>              The cylindrical patch is a 3D geometry that may be used,
        !!              for example, in setting up a cylindrical solid boundary
        !!              confinement, like a blood vessel. The geometry of this
        !!              patch is well-defined when the centroid, the radius and
        !!              the length along the cylinder's axis, parallel to the x-,
        !!              y- or z-coordinate direction, are provided. Please note
        !!              that the cylindrical patch DOES allow for the smoothing
        !!              of its lateral boundary.
        !!      @param patch_id is the patch identifier
    subroutine s_cylinder(patch_id) ! --------------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< Generic loop iterators

        ! Transferring the cylindrical patch's centroid, length, radius,
        ! smoothing patch identity and smoothing coefficient information
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        length_x = patch_icpp(patch_id)%length_x
        length_y = patch_icpp(patch_id)%length_y
        length_z = patch_icpp(patch_id)%length_z
        radius = patch_icpp(patch_id)%radius
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Computing the beginning and the end x-, y- and z-coordinates of
        ! the cylinder based on its centroid and lengths
        x_boundary%beg = x_centroid - 0.5d0*length_x
        x_boundary%end = x_centroid + 0.5d0*length_x
        y_boundary%beg = y_centroid - 0.5d0*length_y
        y_boundary%end = y_centroid + 0.5d0*length_y
        z_boundary%beg = z_centroid - 0.5d0*length_z
        z_boundary%end = z_centroid + 0.5d0*length_z

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smearing of the cylindrical patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the cylinder covers a particular cell in the
        ! domain and verifying whether the current patch has the permission
        ! to write to that cell. If both queries check out, the primitive
        ! variables of the current patch are assigned to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then

                        if (length_x /= dflt_real) then
                            eta = tanh(smooth_coeff/min(dy, dz)* &
                                       (sqrt((cart_y - y_centroid)**2 &
                                             + (cart_z - z_centroid)**2) &
                                        - radius))*(-0.5d0) + 0.5d0
                        elseif (length_y /= dflt_real) then
                            eta = tanh(smooth_coeff/min(dx, dz)* &
                                       (sqrt((x_cc(i) - x_centroid)**2 &
                                             + (cart_z - z_centroid)**2) &
                                        - radius))*(-0.5d0) + 0.5d0
                        else
                            eta = tanh(smooth_coeff/min(dx, dy)* &
                                       (sqrt((x_cc(i) - x_centroid)**2 &
                                             + (cart_y - y_centroid)**2) &
                                        - radius))*(-0.5d0) + 0.5d0
                        end if

                    end if

                    if ((((length_x /= dflt_real .and. &
                           (cart_y - y_centroid)**2 &
                           + (cart_z - z_centroid)**2 <= radius**2 .and. &
                           x_boundary%beg <= x_cc(i) .and. &
                           x_boundary%end >= x_cc(i)) &
                          .or. &
                          (length_y /= dflt_real .and. &
                           (x_cc(i) - x_centroid)**2 &
                           + (cart_z - z_centroid)**2 <= radius**2 .and. &
                           y_boundary%beg <= cart_y .and. &
                           y_boundary%end >= cart_y) &
                          .or. &
                          (length_z /= dflt_real .and. &
                           (x_cc(i) - x_centroid)**2 &
                           + (cart_y - y_centroid)**2 <= radius**2 .and. &
                           z_boundary%beg <= cart_z .and. &
                           z_boundary%end >= cart_z)) &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                    end if

                end do
            end do
        end do

    end subroutine s_cylinder ! --------------------------------------------

    !>      The swept plane patch is a 3D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, at an angle with respect to the axes of the
        !!              Cartesian coordinate system. The geometry of the patch is
        !!              well-defined when its centroid and normal vector, aimed
        !!              in the sweep direction, are provided. Note that the sweep
        !!              plane patch DOES allow the smoothing of its boundary.
        !!      @param patch_id is the patch identifier
    subroutine s_sweep_plane(patch_id) ! -----------------------------------

        integer, intent(IN) :: patch_id

        integer :: i, j, k !< Generic loop iterators

        ! Transferring the centroid information of the plane to be swept
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        smooth_coeff = patch_icpp(patch_id)%smooth_coeff

        ! Obtaining coefficients of the equation describing the sweep plane
        a = patch_icpp(patch_id)%normal(1)
        b = patch_icpp(patch_id)%normal(2)
        c = patch_icpp(patch_id)%normal(3)
        d = -a*x_centroid - b*y_centroid - c*z_centroid

        ! Initializing the pseudo volume fraction value to 1. The value will
        ! be modified as the patch is laid out on the grid, but only in the
        ! case that smearing of the sweep plane patch's boundary is enabled.
        eta = 1d0

        ! Checking whether the region swept by the plane covers a particular
        ! cell in the domain and verifying whether the current patch has the
        ! permission to write to that cell. If both queries check out, the
        ! primitive variables of the current patch are written to this cell.
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    if (grid_geometry == 3) then
                        call s_convert_cylindrical_to_cartesian_coord(y_cc(j), z_cc(k))
                    else
                        cart_y = y_cc(j)
                        cart_z = z_cc(k)
                    end if

                    if (patch_icpp(patch_id)%smoothen) then
                        eta = 5d-1 + 5d-1*tanh(smooth_coeff/min(dx, dy, dz) &
                                               *(a*x_cc(i) + &
                                                 b*cart_y + &
                                                 c*cart_z + d) &
                                               /sqrt(a**2 + b**2 + c**2))
                    end if

                    if ((a*x_cc(i) + b*cart_y + c*cart_z + d >= 0d0 &
                         .and. &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i, j, k))) &
                        .or. &
                        patch_id_fp(i, j, k) == smooth_patch_id) &
                        then

                        call s_assign_patch_primitive_variables(patch_id, i, j, k)

                    end if

                end do
            end do
        end do

    end subroutine s_sweep_plane ! -----------------------------------------


end module m_create_patches
