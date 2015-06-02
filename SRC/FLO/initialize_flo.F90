MODULE initialize_flo
#include "ibm.h"
	use precision 
	use constants 
	!use scalar_init
	use general_funcs
	use randomno
	use hard_sphere
	use postproc_funcs
	use initialize
	use fftw3_interface
	!use mypost_process
	use init_turb, only : specinit, calc_initial_turb, tke_i, eddy_time_i, tau_eta_i, u_prime
	use rand_number

	implicit none 
	!Read input files for initial velocity field and simulation
	!parameters. Check initial flow field for continuity and 
	!correct for any divergence. 

	integer, parameter, private :: nbin_max = 200, mis_gof=120
	integer, private :: gofunit, gofavgunit
	real(prcn), private ::  gofr(nbin_max), rho_est(nbin_max)
	LOGICAL, private :: rescaling
CONTAINS
	subroutine initflo

use mypost_process, only : velocity_output, uiui_correlation

		use dependent_functions
		use restart_funcs
		use global_data
		use geom_init
		!use nlarrays, ONLY  : uf1, ur1
		use nlmainarrays, Only : ubcp, pbcp !, nlbcp, onlbcp , ubc, nlbc, onlbc, pbc
		use maternmod
		use collision_mod

		!^^^ Mohammad 10-16-2009 ^^^
		!use mypost_process
		!use string_funcs
		!use init_turb, only : avr, tke, turb_uf
		!use mypost_process
		!use steptheflow
		!---------------------------

		implicit none  

		!-----------------------------------------------------------------------
		!local variables

		real(prcn) ::  jt(ndim,ndim), flagu,dum, gammaloc, ldomtmp(3), upi_int(ndim), uchar_tmp
		integer :: ii,i,j,k,l,n,m, mytmp, mbox,idim, partstart, partend, iphs
		real(prcn) ::   umax_tmp, lchar, xli(ndim), xlo(ndim), vfrac_mean !,umeanslip
		real(prcn) :: vbox
		real(prcn) :: dt_turb


real(prcn), allocatable :: center(:,:)
real(prcn) :: dist1, dist2

real(prcn), allocatable, dimension(:,:,:) :: f_1,f_2,f_3,f_4
complex(prcn), allocatable, dimension(:,:,:) :: fc_1,fc_2,fc_3,fc_4
		
		if (irestart==0) call init_rand

		first_pass = .true.
		soft_start = .true.
		move_particles = .false.
		!if (imove.eq.1)soft_start = .false.
		base_state_saved = .false.
		compute_auto_corr = .false.
		rhoofs0 = one
		stepbody=1
		gauss_u = .false.
		gauss_p = .false.
		gauss_phi = .false.
		gauss_p = .false.
		if (only_dem) move_particles = .true.
		if (.not.mean_vel_to_particles) movingcv = .false.
		if (imove.eq.1) movingcv = .false.

		call set_interpolation_scheme(2)
		! do it before allocating memory
		! for gstencil and rest of the arrays required for interpolation.
		! This iniliazes order variable required for allocating other
		! arrays. do it for the largest shcmes that is gonna be used. 

		if (irestart.eq.0) then
			count_restart = 0
		endif
    ! In this routine, if MY is undefined, it will be set based on LYBYD and dbydx. if input_type is random, then LYBYD is based on the average diameter and dbydx is based on the smallest diameter.

		call generate_configuration

		if (I_AM_NODE_ZERO) call screen_separator(80,'#')
		BARRIER(comm_cart_2d)

		vbox = doml(1)*doml(2)*doml(3)

		nphsc2 = nphases*(nphases-1)/2
		partstart = 1
		do iphs = 1, nphases
			phase_array(iphs)%pstart = partstart
			if (irestart.eq.0) phase_array(iphs)%dia = radbdy(partstart)*two*dx
			!print*,'radius = ', iphs, phase_array(iphs)%dia, dx
			phase_array(iphs)%volfrac = (phase_array(iphs)%npart)*pi*(phase_array(iphs)%dia**3.d0)/(6.d0*vbox)

			partend = partstart + phase_array(iphs)%npart-1 
			phase_array(iphs)%pend = partend
			do m = PARTSTART,PARTend
				part_array(m)%iphs = iphs
			enddo
			partstart = partend + 1
		enddo
		if (I_AM_NODE_ZERO) then
			if ((trim(input_type).eq.'random').and.(trim(psd_type).eq.'bidisp')) then
				write (*,'(A30,I6, 2x,g17.8)') 'VOLFRAC RATIO FOR BIDISP CASE = ', myid, phase_array(2)%volfrac/phase_array(1)%volfrac
			endif
		endif

		mean_volfrac = zero
		do iphs = 1, nphases
			mean_volfrac = mean_volfrac + phase_array(iphs)%volfrac
		enddo



		if (igeometry==1 .and. from_post) return


		!XLENGTH = doml(1)
		!ZLENGTH = doml(3)
		!YLENGTH = doml(2)

		!if (from_post.and.POST_NO_FLOW_MEM_ALLOC) return

		call initiate_fftw3

		if (.not.allocated(u)) then
			call alloc_mem !allocate memory for main data 
			!ALSO allocate SCALAR 
			!if (ISCALON.EQ.1) call alloc_scalar_mem
		endif

#if !PARALLEL
		do i=local_o_start(1), local_o_end(1)
			ii = i-local_o_start(1)+1

			if (i>local_ni(1)/2+1) then
				write (*,*) myid, "PROBLEM IN THE WX"
				PARALLEL_FINISH()
				stop
			endif

			wx(ii) = dcmplx( zero, twopi*dble(i-1) / (global_n(1)*dx) )
		enddo

		do i=local_o_start(2), local_o_end(2)
			ii = i-local_o_start(2)+1
			if (i<=global_n(2)/2) then
				wy(ii) = dcmplx( zero,  twopi*dble(i-1)             / (global_n(2)*dy) )
			else
				wy(ii) = dcmplx( zero, -twopi*dble(global_n(2)+1-i) / (global_n(2)*dy) )
			endif
		enddo

		do i=local_o_start(3), local_o_end(3)
			ii = i-local_o_start(3)+1
			if (i<=global_n(3)/2) then
				wz(ii) = dcmplx( zero,  twopi*dble(i-1)             / (global_n(3)*dz) )
			else
				wz(ii) = dcmplx( zero, -twopi*dble(global_n(3)+1-i) / (global_n(3)*dz) )
			endif
		enddo

		!	calculate diffusion coefficient at each grid point
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					!new w2(i,j,k)=-dreal(wy(j)*wy(j)+wz(k)*wz(k))
					w2(i,j,k)= -dreal(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k)) !-dreal(wx(i)*wx(i)+wz(k)*wz(k))
				enddo
			enddo
		enddo

#else
		!WAVE NUMBER WY IS ASSOCIATED WITH X-DIRECTION IN REAL SPACE IN A TRANSPOSED FORM
		do i=local_o_start(1), local_o_end(1)
			ii = i-local_o_start(1)+1

			if (i>global_n(1)/2+1) then
				write (*,*) myid, "PROBLEM IN THE WX"
				PARALLEL_FINISH()
				stop
			endif

			wy(ii) = dcmplx( zero, twopi*dble(i-1) / (global_n(1)*dx) )
		enddo

		!WAVE NUMBER WZ IS ASSOCIATED WITH Y-DIRECTION IN REAL SPACE IN A TRANSPOSED FORM
		do i=local_o_start(2), local_o_end(2)
			ii = i-local_o_start(2)+1
			if (i<=global_n(2)/2) then
				wz(ii) = dcmplx( zero,  twopi*dble(i-1)             / (global_n(2)*dy) )
			else
				wz(ii) = dcmplx( zero, -twopi*dble(global_n(2)+1-i) / (global_n(2)*dy) )
			endif
		enddo

		!WAVE NUMBER WX IS ASSOCIATED WITH Z-DIRECTION IN REAL SPACE IN A TRANSPOSED FORM
		do i=local_o_start(3), local_o_end(3)
			ii = i-local_o_start(3)+1
			if (i<=global_n(3)/2) then
				wx(ii) = dcmplx( zero,  twopi*dble(i-1)             / (global_n(3)*dz) )
			else
				wx(ii) = dcmplx( zero, -twopi*dble(global_n(3)+1-i) / (global_n(3)*dz) )
			endif
		enddo

		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					w2(i,j,k)=-dreal(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k))
				enddo
			enddo
		enddo
#endif

		!-----------------------------------------------------------------------
		! calculate phase-shift coefficients
		!if (aliasflag.eq.1) then 
		!	do k=1,mz
		!		do j=1,my
		!			shiftyz(j,k)=EXP((wy(j)*dy+wz(k)*dz)/three)
		!		enddo
		!	enddo
		!endif

		!-----------------------------------------------------------------------


		!-----------------------------------------------------------------------
		!	define maximum wavenumber for spectrum truncation

		!wmax2=w2(my2,1)
		!wy(my2)=czero
		!wz(mz/2+1)=czero

		if (irestart==0) then
			if (iturbon) call specinit
		else
			if (iturbon) call calc_initial_turb
		endif

! call uiui_correlation
!stop

#if 0
#if !PARALLEL
		if (mod(global_n(1),2)==0) wx(global_n(1)/2+1) = czero
		if (mod(global_n(2),2)==0) wy(global_n(2)/2+1) = czero
		if (mod(global_n(3),2)==0) wz(global_n(3)/2+1) = czero

		!	calculate diffusion coefficient at each grid point
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					!new w2(i,j,k)=-dreal(wy(j)*wy(j)+wz(k)*wz(k))
					w2(i,j,k)= -dreal(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k)) !-dreal(wx(i)*wx(i)+wz(k)*wz(k))
				enddo
			enddo
		enddo
#else
		if (mod(global_n(1),2)==0) then
			i = global_n(1)/2+1
			if ( local_o_start(1) <= i .and. i <= local_o_end(1) ) then
				ii = i-local_o_start(1)+1
				wy(ii) = czero
			endif
		endif

		if (mod(global_n(2),2)==0) then
			i = global_n(2)/2+1
			if ( local_o_start(2) <= i .and. i <= local_o_end(2) ) then
				ii = i-local_o_start(2)+1
				wz(ii) = czero
			endif
		endif

		if (mod(global_n(3),2)==0) then
			i = global_n(3)/2+1
			if ( local_o_start(3) <= i .and. i <= local_o_end(3) ) then
				ii = i-local_o_start(3)+1
				wx(ii) = czero
			endif
		endif

		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					w2(i,j,k)=-dreal(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k))
				enddo
			enddo
		enddo
#endif
#endif

		call grid_nodes_insphere

		if (I_AM_NODE_ZERO) then
			write (*,'(A)') 'OUT OF GRID_NODES_INSPHERE'
			write (*,'(2(2x,A25,i10))') 'COUNT_SOLID = ', count_solid, 'COUNT_FLUID = ', count_fluid
			write (*,'(A10,2(2x,g17.8))') 'MAXVOLFRAC = ', maxvolfrac, count_solid
		endif
		
		grav(:) = zero
		if (impose_grav) then
			set_mpg = .true.

			grav(1) = archno * vis**2 / char_length**3 * rhof / abs(rhos-rhof)
			grav(2) = zero
			grav(3) = zero

			mpg(:) = -grav(:) * rhof * (1+maxvolfrac*(rhos/rhof-1))

			!mpg(:) = (-0.2638240D-02) * cos(pi*flo_ang(1:3)/180.d0)
			!grav(:) = - mpg(:) / (rhof*(one-maxvolfrac) + rhos*maxvolfrac)
			!grav(1:3) = -0.2638240D-02*cos(pi*flo_ang(1:3)/180.d0)
			!mpg(:) = (rhof*(one-maxvolfrac) + rhos*maxvolfrac)*grav(:)
			!vis = DSQRT(ABS(rhos/rhof-one)* sqrt(dot_product(grav(:),grav(:))) * char_length**3.d0/archno)
			!if (I_AM_NODE_ZERO) write (*,'(2(A25,2x,g17.8))') 'TO ATTAIN AR. NUMBNER = ', archno, 'CHANGING VISCOSITY TO = ', vis

			uchar_tmp = sqrt( sqrt(dot_product(grav(:),grav(:))) *(rhos/rhof-1)*char_length) * (1-maxvolfrac)
			uchar(1:3) = uchar_tmp*cos(pi*flo_ang(1:3)/180.d0)
			ucharmod = SQRT(uchar(1)**2.d0 + uchar(2)**2.d0+uchar(3)**2.d0)    
			fsslip(1:ndim) = uchar(1:ndim)/(one-maxvolfrac)

			fsslipmod = sqrt(dot_product(fsslip,fsslip))

			if (I_AM_NODE_ZERO) then
				write (*,"(1a,3d15.7)") "Ar NUMBER                  = ", archno
				write (*,"(1a,3d15.7)") "GRAVITATIONAL ACCELERATION = ", grav(:)
				write (*,"(1a,3d15.7)") "MPG                        = ", mpg(:)

			endif

		elseif (zero_slip) then
			set_mpg = .true.

			mpg(:) = zero
			grav(:) = zero

			if (iturbon) then
				uchar_tmp = u_prime
				uchar(1:3) = uchar_tmp 
				ucharmod = uchar_tmp
				fsslip(1:ndim) = zero
				fsslipmod = sqrt(dot_product(fsslip,fsslip))


				if (I_AM_NODE_ZERO) then
					write (*,*) "SCALES ARE SELECTED BASED ON INITIAL ISOTROPIC TURBULENCE"
					write (*,"(1a,1d15.7)") "VELOCITY SCALE = ", u_prime
					write (*,"(1a,1d15.7)") "TIME SCALE     = ", eddy_time_i
				endif
			elseif (Ret>small_number) then
				gran_temp = (ReT*vis/char_length)**2.d0
				grant_i = gran_temp

				uchar_tmp  = sqrt(gran_temp)
				uchar(1:3) = uchar_tmp 
				ucharmod = uchar_tmp !SQRT(uchar(1)**2.d0 + uchar(2)**2.d0+uchar(3)**2.d0)    
				fsslip(1:ndim) = zero !uchar(1:ndim)/(one-maxvolfrac)
				fsslipmod = sqrt(dot_product(fsslip,fsslip))

				if (I_AM_NODE_ZERO) then
					write (*,*) "SCALES ARE SELECTED BASED ON INITIAL GRANULAR TEMPERETURE"
					write (*,"(1a,1d15.7)") "VELOCITY SCALE = ", ucharmod
					write (*,"(1a,1d15.7)") "TIME SCALE     = ", char_length/ucharmod
				endif
			endif
		else
			!uchar_tmp = ((Re*vis)/char_length)!*cos(pi*flo_ang(1:3)/180.d0)
			! New formulation to make IBM GI

			fsslipmod = one
			uchar_tmp = (one-maxvolfrac)*fsslipmod
			vis = uchar_tmp*char_length/Re 
			if (I_AM_NODE_ZERO) write (*,*)'CHANGING VISCOSTY TO GET UCHAR = 1. NEW VIS = ', vis

			uchar_tmp = ((Re*vis)/char_length)
			uchar(:) = uchar_tmp*cos(pi*flo_ang(:)/180.d0)
			ucharmod = sqrt(dot_product(uchar(:),uchar(:)))
			fsslip(:) = uchar(:)/(one-maxvolfrac)
			!dchar = dia_phys
		endif

		if (ReT.gt.SMALL_NUMBER) gran_temp = (ReT*vis/char_length)**2.d0

 		! Assign mean and fluctuating velocities to particles
		if (irestart.eq.0) then
			if (trim(input_type).eq."random") then
				call generate_particle_velocities
			elseif (trim(input_type).eq."default") then
				call generate_particle_velocities
			elseif (trim(input_type).eq."lubtest") then
				if (I_AM_NODE_ZERO) then
					write (*,*) 'LUBRICATION TEST CASE CASE. CHARACTERISTIC VEL SCALE IN RE IS REL. VELOCITY.'
					write (*,*) 'OPPOSITE VELOCIIES ARE ASSIGNED.'
				endif
				velbdy(1,1:ndim) = uchar(1:ndim)/two
				velbdy(2,1:ndim) = -uchar(1:ndim)/two
			else
				if (mean_vel_to_particles) then
					do idim=1,ndim
						do m = 1, nbody
							velbdy(m,idim) =  -uchar(idim)/(one-maxvolfrac)
						enddo
					enddo
				else
					do m = 1, nbody
						velbdy(m,1:ndim) = zero
					enddo
				endif
			endif
    
			do idim = 1, ndim
				if (nbody.gt.0) then
					usmean_des(idim) = SUM(velbdy(1:nbody, idim))/real(nbody, prcn)
				else
					usmean_des(idim) = zero
				endif
			enddo
    
			! Assign mean velocity to the fluid.             
			if (mean_vel_to_particles) then
				ufmean_des(:) = zero
			elseif (zero_slip) then
				ufmean_des(:) = zero
			else
				ufmean_des(:) = uchar(:)/(one-maxvolfrac)
			endif

			if (trim(input_type).eq.'lubtest') ufmean_des(1:ndim) = zero
    
			! Transform velocities to get homogeneous boundary conditions on particle surface
			do idim = 1, ndim
				!!$do m = 1, nbody
				!!$velbdy(m,idim) = (velbdy(m,idim) - solid_vel(idim))/fsslipmod
				!!$enddo
				!!$usmean_des(idim) = (usmean_des(idim)-solid_vel(idim))/fsslipmod
				!!$ufmean_des(idim) = (ufmean_des(idim)-solid_vel(idim))/fsslipmod
          
				! Compute the volumetric mean velocity
				if (.not.set_mpg) then
					umean(idim) = (one-maxvolfrac)*ufmean_des(idim) + maxvolfrac*usmean_des(idim)
				else
					!umean(idim) = zero
					umean(idim) = (one-maxvolfrac)*ufmean_des(idim) + maxvolfrac*usmean_des(idim)
				endif
			enddo
    
			frame_accln(1:ndim) = zero
			frame_vel(1:ndim) = zero
			frame_pos(1:ndim) = one
			call INPUT_CHECK
		endif


#if 0
		do i=1, 1
			do j=1, 1l
				do k=1, 1
					do m=1, nbody
						write (1,"(4d15.7,1i)") xc(m,1)/dbydx+(i-1)*lybyd, &
					&									xc(m,2)/dbydx+(j-1)*lybyd, &
					&									xc(m,3)/dbydx+(k-1)*lybyd, radbdy(m)/dbydx, 1
					enddo
				enddo
			enddo
		enddo

		write (1,*) "zone"
		do i=1, 1
			do j=1, 1
				do k=1, 1
					do m=1, nbody
						if (xc(m,2)/dbydx<=lybyd/4.or.xc(m,2)/dbydx>=lybyd*3/4) &
					&	write (1,"(4d15.7,1i)") xc(m,1)/dbydx+(i-1)*lybyd, &
					&									xc(m,2)/dbydx+(j-1)*lybyd, &
					&									xc(m,3)/dbydx+(k-1)*lybyd, radbdy(m)/dbydx, 1
					enddo
				enddo
			enddo
		enddo

		close(1)
		PARALLEL_FINISH()
		stop
#endif

#if 0
		if (I_AM_NODE_ZERO) then
			write (*,*) 'WRITING IMMERSED BODY DATA'
			open (unit=1,file=trim(run_name)//'_sphr_center.inp',status='replace', Action='write')
			!write (1,"(1a,1d15.7)") "VOLFRAC=", maxvolfrac
			!write (1,"(1i,1a)") nbody, " atoms"
			!write (1,"(1i,1a)") nphases, " atom types"
			!write (1,"(1a)") "atom_style = sphere"
			!write (1,"(2f8.2,1A)") zero, lybyd, " xlo xhi"
			!write (1,"(2f8.2,1A)") zero, lybyd, " ylo yhi"
			!write (1,"(2f8.2,1A)") zero, lybyd, " zlo zhi"
			!write (1,"(1a)") "atoms"
			!write (1,"(1a)") " "
			do i=1, nbody
				!write (1,"(2i,5d15.7)") i, 1, dia_phys, 1000.0, xc(i,:)/dbydx
				write (1,"(4d15.7)") xc(i,:)/dbydx, radbdy(i)/dbydx
			enddo
			!write (1,"(1a)") "Velocities"
			!write (1,"(1a)") " "
			!do i=1, nbody
			!	write (1,"(1i,6d15.7)") i, velbdy(i,:), 0d0, 0d0, 0d0
			!enddo
			close(1)
		endif
		PARALLEL_FINISH()
		stop
#endif


		dtorig = dt
		urf = 1.0d0 ! momentum equation under-relaxation factor
		prf = 1.0d0!0.8d0
		!r = 0.5*dia_phys/doml(2)*my
		upi = zero 

#if 0
		if (I_AM_NODE_ZERO.AND.nbody.gt.0) then
			open (unit=2001,file=trim(run_name)//'_sphr_vel_out.dat',form='formatted',status='unknown')
			write (2001,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ', ' "UX" '
			do m=1,nbody
				write (2001,'(6(2x,f12.8))') xc(m,1), xc(m,2), xc(m,3), velbdy(m,1), velbdy(m,2), velbdy(m,3)
			enddo
			close(2001, status='keep')
		endif
#endif

		if (nphases.gt.0) then
			r = 0.5*phase_array(nphases)%dia/doml(2)*my ! Resolve the largest sphere and store in the array of first phase
			do iphs = 1, 1!nphases
				!r = half*(phase_array(iphs)%dia)/doml(2)*my
				if (I_AM_NODE_ZERO) write (*,*)'RADIUS IN QUADGENER = ', r
				nullify(phase_array(iphs)%bndpts)
				if (allocated(xs)) deallocate(xs)
				if (allocated(cd)) deallocate(cd)
				call quadgener(r,dr,f1,f2)

				nrpr = nbnd
				if (I_AM_NODE_ZERO) then
					write (*,'(A,2(2x,I5))') 'NO. OF BOUNDARY OPINTS FOR PHASE ', iphs, nbnd           
					write (ounit,'(A,2(2x,I5))') 'NO. OF BOUNDARY OPINTS FOR PHASE ', iphs, nbnd           
				endif
 
				if (associated(phase_array(iphs)%bndpts)) deallocate(phase_array(iphs)%bndpts)
				if (nbnd.GT.0) allocate(phase_array(iphs)%bndpts(ndim,nbnd))
    
				phase_array(iphs)%nbnd = nbnd
				phase_array(iphs)%nrpr = nrpr

				do i = 1, nbnd
					do idim = 1, ndim
						phase_array(iphs)%bndpts(idim,i) = xs(idim,i)
					enddo
				enddo
			enddo
			if (allocated(xs)) deallocate(xs)
		endif

		do m = 1, nbody
			radibdy(m) = radbdy(m)-dr
			radobdy(m) = radbdy(m)+dr
			rado2bdy(m) = radbdy(m)+two*dr
		enddo

		do m= 1, nbody
			iphs = 1 !part_array(m)%iphs
			nrpr = phase_array(iphs)%nrpr
			nullify(part_array(m)%if_rev)
			nullify(part_array(m)%if_drag)

			if (associated(part_array(m)%if_rev))  deallocate(part_array(m)%if_rev)
			if (associated(part_array(m)%if_drag)) deallocate(part_array(m)%if_rev)

			if (nrpr.GT.0) then
				allocate(part_array(m)%if_rev(nrpr))
				allocate(part_array(m)%if_drag(nrpr))
			endif
		enddo

		if (irestart.EQ.1) then
			call read_restart
			call INPUT_CHECK
		endif
 		voldom = doml(1)*doml(2)*doml(3)

		if (zero_slip) then
			if (iturbon) then
				t_conv = eddy_time_i
				tendused = 4*t_conv
			elseif (Ret>small_number) then
				t_conv = char_length/ucharmod
				tendused = 4*t_conv
			endif
		else
			uchar_tmp = SQRT(uchar(1)**2.d0 + uchar(2)**2.d0+uchar(3)**2.d0)/(one-maxvolfrac)
			umeanslip = uchar_tmp

			t_conv = (char_length)/(uchar_tmp)
			tendused = tend*doml(1)/(uchar_tmp)
		endif

		if (I_AM_NODE_ZERO) then
			write (*,"(1a,10d15.7)") 'CHAR LENGTH, UCHAR_TEMP = ', char_length, uchar_tmp
			write (*,"(1a,10d15.7)") 'UCHAR(1:3), UCHARMOD    = ', uchar, ucharmod
			write (*,"(1a,10d15.7)") 'FSSLIP(1:3), FSSLIPMOD  = ', fsslip(1:ndim), fsslipmod
			write (*,"(1a,3d15.7)") "GRAVITY = ", grav(:)
			write (*,"(1a,3d15.7)") "MPG     = ", mpg(:)
			if (zero_slip.and.Ret>small_number) then
				write (*,"(1a,2d15.7)") "ReT, T/W^2= ", ReT, gran_temp/ ucharmod**2
			endif
			write (*,'(A,2(2x,g12.6))') 'UCHAR AND MAXVOLFRAC = ', UMEANSLIP, maxvolfrac
		endif

		t_min = t_conv
		umax_tmp = MAX(uchar_tmp, MAXVAL(ABS(velbdy(1:nbody,1:ndim))))

		if (I_AM_NODE_ZERO) write (*,'(A25,g12.6)') 'MAXIMUM VELOCITY  = ', umax_tmp

		dt_tmp_conv = (dx)/(umax_tmp)!/100
		dt_tmp = dt_tmp_conv 

		!lchar = dchar!*(one-maxvolfrac)
		lchar  = char_length
		t_vis =(lchar*lchar)/vis

		diff_ts = .false.

		dt_tmp_vis = large_number
		dt_tmp_diff = large_number

		if (t_vis-t_min.LT.small_number) then 
			if (I_AM_NODE_ZERO) print*,'T_VIS-T_MIN = ', T_VIS-T_MIN
			t_min = t_vis
			diff_ts = .true.
			lchar = ((one-maxvolfrac)*doml(1)*doml(2)*doml(3))**(one/three)
			!lchar = lchar*doml
			tendused = 0.2*tend*lchar*lchar/vis
			dt_tmp = dt_tmp_vis
		endif
		lchar = dx 

		if (Re.lt.one) then
			!dt_tmp_vis = (lchar*lchar*(one-maxvolfrac))/vis
			!dt_tmp_vis = t_vis/50.d0 !(lchar*lchar*(one-maxvolfrac))/vis
			!dt_tmp_vis = Re*dbydx*dt_tmp_conv/(cfl*vfl) 
			!vfl is the number
			! of time steps used to resolve the viscous time scale
		endif

		if (I_AM_NODE_ZERO) write (*,'(A25,g12.6)') 'RE = tvis/tcon = ', t_vis/t_conv

		!t_grav = DSQRT(two*dchar/9.81)
		t_grav = DSQRT(two*char_length/sqrt(dot_product(grav,grav)))
		dt_tmp_grav = LARGE_NUMBER
		if (impose_grav) dt_tmp_grav = DSQRT(two*dx/sqrt(dot_product(grav,grav)))

		!dt_tmp_vis = (vfl*cfl*lchar*lchar)/vis
		!t_vis =( 4.d0*r*r)/vis
		!endif

		dt_tmp_diff = LARGE_NUMBER
		t_diff = LARGE_NUMBER 
		if (iturbon) dt_turb = dx * mx/dble(mx_iso) /20/sqrt(tke_i)

#if 0
		if (iscalon.EQ.1) then
			call INITSCAL ! To intialize the Scalar Field
			!dt_tmp_diff = dt_tmp_diff*(1.d0-maxvolfrac) 
		else
			dt_tmp_diff = LARGE_NUMBER
			t_diff = LARGE_NUMBER 
		endif
#endif

		dtcoll = LARGE_NUMBER
		if (nbody.gt.0) then
			XLENGTH = doml(1)
			ZLENGTH = doml(3)
			YLENGTH = doml(2)
			deS_EN_INPUT(:) = coeff_rest
			deS_ET_INPUT(:) = zero

			call des_time_march(.true.)
		endif

		do m=1,nbody            ! loop over bodies
			call update_nrpr_array(m)
		enddo

		if (debug_check) then
			if (I_AM_NODE_ZERO) then
				call screen_separator(40,'^')				
				write (*,*) "PARTICLE ID, DRAG_ACTIVE, NRPR_ACTIVE"
				do m=1, nbody
					write (*,*) m, part_array(m)%drag_active, part_array(m)%nrpr_active
				enddo
				call screen_separator(40,'-')
			endif			
		endif

		if (irestart.eq.0) then
			!if (imove==1 .or. move_particles) then
			!	dtcoll = dtsolid_orig * 50 !new
			!endif
			S_TIME = ZERO 
			AUTO_CORR_SEP_TIME = ZERO
			if (I_AM_NODE_ZERO) then
				print*, 'MAXIMUM OVERLAP = ', OVERLAP_MAX
				!write (*,'(A25,2x,i5)')'default # of Rev pts',nrpr
			endif
    
			if (I_AM_NODE_ZERO) then
				write (*,'(A25,2x,g12.5)')'DT_CONV = ', DT_TMP_CONV
				write (*,'(A25,2x,g12.5)')'DT_VIS = ', DT_TMP_VIS
				write (*,'(A25,2x,g12.5)')'DT_GRAV = ', DT_TMP_grav
				write (*,'(A25,2x,g12.5)')'DT_DIFF = ', DT_TMP_DifF
				if (iturbon) write (*,'(A25,2x,2g12.5)')'DT_TURB = ', dx/20/cfl/sqrt(tke_i)
				write (*,'(A25,2x,g12.5)')'DT_COLLISIONAL = ', DTCOLL
			endif
			dt = MIN(DT_TMP_CONV, DT_TMP_VIS, DT_TMP_grav,DT_TMP_DifF)

			dt = cfl*dt
			if (iturbon) dt = MIN(dt, dt_turb)
			dt = MIN(dt, DTCOLL)
    
			if (I_AM_NODE_ZERO) write (*,'(A25,2x,g12.5)')'DT CHOSEN = ', DT

			cf = -1.d0/dt
			cforig = cf
			if (I_AM_NODE_ZERO) write (*,*) 'cf=',cf
			if (I_AM_NODE_ZERO) write (ounit,'(2x,A,g12.5)') 'cf=',cf

			ferror = one
			fold = zero
    
			do iphs = 1, nphases
				phase_array(iphs)%ferror = one
				phase_array(iphs)%fold = zero
				phase_array(iphs)%ferror_array(1:nerr_steps) = one
				phase_array(iphs)%ferror_hist = one
				if (imove.eq.1) then
					phase_array(iphs)%grant_error = one
					phase_array(iphs)%grant_old = gran_temp
					phase_array(iphs)%grant_array(1:nerr_steps) = one
					phase_array(iphs)%gran_error_hist = one
				endif
			enddo

			ferror_array(:) = one
			ferror_hist = one
			source_hydro(:) = zero 
		endif
 
		!new ramp_frac_time = ramp_frac_steps*ucharmod*dt/(float(mx-1)*dx)
		ramp_frac_time = ramp_frac_steps*ucharmod*dt/(float(global_n(1))*dx)
		if (I_AM_NODE_ZERO) write (*,'(A25,2x,g12.5)')'ramping time steps ', ramp_frac_steps
		if (I_AM_NODE_ZERO) write (*,'(A25,2x,g12.5)')'ramp_frac_time ', ramp_frac_time

		mesh_vel(1:ndim) = usmean_des(1:ndim)
		!if (.NOT.XPERIODIC) then 
		!	do i=1,mx-mbuffer
		!		g(i)=one
		!	enddo
		!	do i=mx-mbuffer,mx
		!		if (XPERIODIC) print*,'if I AM printING FOR PERIODIC CASE, then COME TO ME IN g(x) CALCULATION' 
		!		j=mx-i
		!		g(i)=(one-dcos(twopi*j/(mbuffer*2.d0)))/two
		!		!g(i)=one
		!	enddo
		!	g(mx)=zero
		!else 
		!	!new g(1:nx+1) = one
		!	g(:) = one
		!endif
		!-----------------------------------------------------------------------
		dx2=dx*dx
		!-----------------------------------------------------------------------
		!	define wavenumbers: introduced by Shankar Subramaniam, 2001

		!write (*,*)'Using new definition of wavenumbers'
		!write (*,*)'(consistent with FFTW)'

		!new do i=1,my2
		!new 	wy(i)=dcmplx(zero,twopi*dble(i-1)/(dble(my)*dy))
		!new enddo



111	FORMAT(e12.4)
112	FORMAT(e12.4,5x,e12.4)

		!-----------------------------------------------------------------------
		!initialize boundary conditions in Fourier space

!!$    do n=1,ndim
!!$       do k=1,mz
!!$          do j=1,my2
!!$
!!$             uin(j,k,n)=u(1,j,k,n)
!!$             uout(j,k,n)=u(mx,j,k,n)
!!$
!!$          enddo
!!$       enddo
!!$    enddo

		if (rk3) then
			!FOR RK3 TIME STEPPING
			coef(1,1)=4.d0/15.d0
			coef(1,2)=4.d0/15.d0 !For Alpha, Diff term
			coef(1,3)=8.d0/15.d0 !for Gamma, NL term at n
			coef(1,4)=zero       !For Zeta, ONL TERM

			coef(2,1)=1.d0/15.d0
			coef(2,2)=1.d0/15.d0
			coef(2,3)=5.d0/12.d0
			coef(2,4)=-17.d0/60.d0

			coef(3,1)=1.d0/6.d0
			coef(3,2)=1.d0/6.d0
			coef(3,3)=3.d0/4.d0
			coef(3,4)=-5.d0/12.d0
			itrmax = 3
		else
			!FOR EULER TIME STEPPING 
			coef(1,1)=half
			coef(1,2)=half
			coef(1,3)=three/two
			coef(1,4)=-half
			itrmax = 1
		endif

#if PARALLEL
		vel_converted = .false.
#endif
		if (I_AM_NODE_ZERO) then
			write (*,'(1a18,1i)') '# of I grid pts = ', global_n(1)
			write (*,'(1a18,1i)') '# of J grid pts = ', global_n(2)
			write (*,'(1a18,1i)') '# of k grid pts = ', global_n(3)
			write (*,'(1a5,1d15.7)') 'dx = ', dx
			write (*,'(1a5,1d15.7)') 'dy = ', dy
			write (*,'(1a5,1d15.7)') 'dz = ', dz
			write (*,'(1a5,1d15.7)') 'dr = ', dr
			write (*,'(1a24,1d15.7)') 'Radius(#of grid pts.) = ', R
			write (*,'(1a24,1d15.7)') 'Re =                    ', Re
			write (*,'(1a24,1i)')     'reversl points (yes=1)  ', dorpr
			write (*,'(1a5,3d15.7)') 'upi =      ', uchar(1:ndim)
			write (*,'(1a5,1d15.7)') 'dia_phys = ', dia_phys
    
			write (*,*) 'Tend = ', tendused
			write (*,'(A,i10)') 'MINIMUM NO. STEPS BASED ON CURRENT TIME STEP =  ', nint(tendused/dt)
			write (*,'(6(2x,A25,g12.5,/))') &
				'Convective Time scale =', t_conv, &
				'Viscous Time scale =', t_vis, &
				'gravity Time Scale =', t_grav,&
				'Diffusive Time scale =', t_diff, &
				'Min. Time scale =', t_min, &
				'Time steP Size =', dt
			!   call write_input(ounit)

			write (ounit,'(3(2x,A25,i4,/),5(2x,A25, g12.5,/),3(2x,A25, i4,/),2x,A25,3(g12.5))') &
				'# of I grid pts = ', global_n(1), &
				'# of J grid pts = ', global_n(2), &
				'# of K grid pts = ', global_n(3), &
				'dx = ', dx, &
				'dy = ', dy, &
				'dz = ', dz, &
				'Radius(#of grid pts.) = ', R, &
				'Re = ', Re, &
				'reversl points (yes=1)', dorpr, &
				'upi = ', uchar(1:3)
			write (ounit,'(5(2x,A25,g12.5,/))') &
				'Convective Time scale =', t_conv, &
				'Viscous Time scale =', t_vis, &
				'Diffusive Time scale =', t_diff, &
				'Min. Time scale =', t_min, &
				'Time step size =', dt
		endif

		if (trim(input_type).eq."single-phase") then
			frmean(1:ndim) = zero
			mpg(1:ndim) = zero
		endif
		!saveitns = int(saveitns/100.d0*dia_phys/(ucharmod*dt))
		!if (I_AM_NODE_ZERO)write (*,*)' NEW SAVE ITNS = ', saveitns


#if 0
!THIS SECTION IS FOR THE 1D AND 2D FOURIER SPECTRUM
write (*,*) "ASSIGNING SINE FUNCTION TO VELOCITY FIELD..."

allocate(f_1(local_ni(1), local_ni(2), local_ni(3)))
allocate(f_2(local_ni(1), local_ni(2), local_ni(3)))
allocate(f_3(local_ni(1), local_ni(2), local_ni(3)))
allocate(f_4(local_ni(1), local_ni(2), local_ni(3)))

allocate(fc_1(local_no(1), local_no(2), local_no(3)))
allocate(fc_2(local_no(1), local_no(2), local_no(3)))
allocate(fc_3(local_no(1), local_no(2), local_no(3)))
allocate(fc_4(local_no(1), local_no(2), local_no(3)))


f_1 = zero
f_2 = zero
f_3 = zero
f_4 = zero

fc_1 = czero
fc_2 = czero
fc_3 = czero
fc_4 = czero

allocate(center(2,2))

 center(1,1) = 9
 center(1,2) = 9

 center(2,1) = 27
 center(2,2) = 27

do j=1, local_ni(2)
	do i=1, local_ni(1)
		f_1(i,j,:) = 0.5*(sin(i*dx)+sin(j*dy))

		dist1 = sqrt( (i-center(1,1))**2 + (j-center(1,2))**2 )
		dist2 = sqrt( (i-center(2,1))**2 + (j-center(2,2))**2 )

		f_2(i,j,:) = f_1(i,j,:)
		f_3(i,j,:) = zero

		if (dist1<=3) then
			f_2(i,j,:) = zero
			f_3(i,j,:) = f_1(9,5,:)
		endif

		if (dist2<=3) then
			f_2(i,j,:) = zero
			f_3(i,j,:) = f_1(27,23,:)
		endif

		f_4(i,j,:) = f_2(i,j,:)+f_3(i,j,:)
	enddo
enddo

open (unit=1, file=trim(run_name)//"_velocity.dat", status="replace")
write (1,*) "zone i=", local_ni(1), " j=", local_ni(2), " f=point"
do j=1, local_ni(2)
	do i=1, local_ni(1)
		write (1,"(10d15.7)") i*dx, j*dy, f_1(i,j,1)
	enddo
enddo

write (1,*) "zone i=", local_ni(1), " j=", local_ni(2), " f=point"
do j=1, local_ni(2)
	do i=1, local_ni(1)
		write (1,"(10d15.7)") i*dx, j*dy, f_2(i,j,1)
	enddo
enddo

write (1,*) "zone i=", local_ni(1), " j=", local_ni(2), " f=point"
do j=1, local_ni(2)
	do i=1, local_ni(1)
		write (1,"(10d15.7)") i*dx, j*dy, f_3(i,j,1)
	enddo
enddo

write (1,*) "zone i=", local_ni(1), " j=", local_ni(2), " f=point"
do j=1, local_ni(2)
	do i=1, local_ni(1)
		write (1,"(10d15.7)") i*dx, j*dy, f_4(i,j,1)
	enddo
enddo

!fourier of f1
 call fftwr2c(f_1(1:local_ni(1),1:local_ni(2),1:local_ni(3)), fc_1(:,:,:))

!fourier of f2
 call fftwr2c(f_2(1:local_ni(1),1:local_ni(2),1:local_ni(3)), fc_2(:,:,:))

!fourier of f3
 call fftwr2c(f_3(1:local_ni(1),1:local_ni(2),1:local_ni(3)), fc_3(:,:,:))

!fourier of f4
 call fftwr2c(f_4(1:local_ni(1),1:local_ni(2),1:local_ni(3)), fc_4(:,:,:))

!compute f4-f3 in fourier
do k=1, local_no(3)
	do j=1, local_no(2)
		do i=1, local_no(1)
			u(i,j,k,1) = fc_4(i,j,k)-fc_3(i,j,k)
		enddo
	enddo
enddo

 call fftwc2r(u(1:local_no(1),1:local_no(2),1:local_no(3),1),  f_2(1:local_ni(1),1:local_ni(2),1:local_ni(3)))

write (1,*) "zone i=", local_ni(1), " j=", local_ni(2), " f=point"
do j=1, local_ni(2)
	do i=1, local_ni(1)
		write (1,"(10d15.7)") i*dx, j*dy, f_2(i,j,1)
	enddo
enddo
 close(1)

write (*,*) global_n(1), global_n(2), global_n(3)
!call velocity_output
stop
#endif

	end subroutine initflo


	subroutine generate_configuration
		use global_data
		use postproc_funcs
		use hard_sphere
		use maternmod
		use collision_mod
		use randomno
		use dependent_functions
		implicit none

		real(prcn) :: rad, temp, dia_gcg1,dia_gcg2,phiavg_gcg, vol_frac2_gcg, vol_tmp

		real(prcn), allocatable, dimension(:,:,:,:) ::  velr
		real(prcn) :: ldomtmp(3), x1,x2, umf0(ndim), rsf(ndim,ndim)

		integer :: i,j,k,n,m,  mytmp, mbox, nbdymax, idim, tmp_unit, iphs, nsim
		real(prcn) :: lambda_p,  rad_tmp, final_vol_frac, pvel_var1, pvel_var2

		real(prcn), dimension(:,:), allocatable :: xc_temp, vtemp
		character*80 :: TEMPCHAR
		character*80 :: filename, filename1

		real(prcn) :: part_sep, dummyx, dummyy, dummyz, dummyr
		integer :: dummyi

		integer, dimension(:), allocatable :: npart

		real(prcn) :: conf, confint, int_dist_avg, int_dist_var, int_dist_sd, int_dist2

		logical, allocatable :: contact(:,:)
		real(prcn) :: max_overlap

		mytmp = my
		mbox = my

		!nbins = 200 
		rescaling = .true.

#if 0
		if (I_AM_NODE_ZERO) then
			gofunit = getnewunit(minunitno,maxunitno)
			open (gofunit,file=trim(run_name)//'_gof_used.dat', form = 'formatted')
			gofavgunit = getnewunit(minunitno,maxunitno)
			open (gofavgunit,file= trim(run_name)//'_gof_avg.dat', form = 'formatted')
		endif
#endif


		if (nbins.GT.nbin_max) then
			write (*,*) "INCREASE nbin_max in initialize flo "
			stop
		endif
    
		!radbdy(1:nbody) = 0.0625
		!allocate the radbdy array to the size of the number of bodies
		if (irestart.eq.0) then 
			if (input_type.eq."mat") then 
				call screen_separator(80,'I')
				write (*,*) "IN MATERN CARD-CORE PROCESS"

				nsim = mis_mat
				if (my.eq.undefined_I) then
					my = lybyd * dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					!new if (my/2.ne.0)my = my + 1
					mz = my
					mx = my !+1
					!mxf = mx 
					!my2 =  my/2+1
					!mx1 = mx !- 1 
					!mx2 = mx1/2+1		
					if (I_AM_NODE_ZERO) print*,'MY UNdeFINED: SETTING TO:', MY
				endif

				if (.not.allocated(gofr_avg)) then
					allocate(gofr_avg(nbins), gofr_mis(nsim, nbins))
					allocate(rad_bin(nbins))
					allocate(int_dist(nsim))

					gofr_avg = zero
					gofr_mis = zero
					rad_bin  = zero
					int_dist = zero
				endif

				doml(2) = lybyd*dia_phys
				dx = doml(2)/(my)
				dy = dx
				dz = dx 
				XLENGTH = doml(1)
				ZLENGTH = doml(3)
				YLENGTH = doml(2)

				!new doml(1) = (mx-1)*dx
				doml(1) = mx*dx
				doml(3) = mz*dz
          
2002			continue
				rad_tmp = half*(real(my))/lybyd
				!rad_tmp = half*(one)/lybyd
				vol_tmp =  phiavg
				lambda_p = -one/(fourthirdpi*(hhat*rad_tmp)**3)

				write (*,"(1A,1D15.7)") "RAD_TMP = ", rad_tmp
				write (*,"(1A,1D15.7)") "LAMBDA_P = ", lambda_p

				if (one.gt.(hhat**3.0)*phiavg) then 
					!lambda_p = lambda_p * log (one-fourthirdpi*(hhat*rad_tmp*2)**3*phiavg)
					lambda_p = lambda_p * log (one-(hhat)**3*phiavg)
					write (*,"(1A,1D15.7)") "LAMBDA_P = ", lambda_p
				else
					write (*,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'
					write (ounit,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'

					write (*,"(1A,1D15.7)") "HHAT = ", hhat
					write (*,"(1A,1D15.7)") "phiavg = ", hhat
					write (*,"(1A,1D15.7)") "(HHAT**3.0)*phiavg = ", (hhat**3.0)*phiavg
					write (*,"(1A,1D15.7)") "LAMBDA_P = ", lambda_p

					PARALLEL_FINISH()
					stop
				endif

				if (I_AM_NODE_ZERO) then
					!new ldomtmp(1) = mx1 !-1
					ldomtmp(1) = mx
					ldomtmp(2:3) = my

					nbdymax = 2*Nint(lambda_p*ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
					write (*,*)'MAX NUMBER OF BODIES IN MATERN : ', nbdymax

					if (allocated(xc_temp)) deallocate(xc_temp)
					allocate(xc_temp(nbdymax,3))

					call matern(3, ibordel, Ldomtmp, lambda_p, rad_tmp, hhat, nsim, nbdymax, xc_temp, nbody)!0.006 
				endif
				BROADCAST_INT(nbody,1,node_zero,comm_cart_2d)
				call alloc_bnd_related
				if (I_AM_NODE_ZERO) then
					do n = 1, ndim
						xc(1:nbody,n) = xc(1:nbody,n) + one
					enddo
					phiavg = vol_tmp
					print*,'NUMBER OF BODIES AFTER MATERN = ', nbody
					write (ounit,*) 'NUMBER OF BODIES AFTER MATERN = ', nbody
					!stop

					xc(1:nbody,1:ndim)  = xc_temp(1:nbody, 1:ndim)
					radbdy(1:nbody)  = rad_tmp

#if 0
					open (unit=2000,file=trim(run_name)//'_sphr_center_mat.dat',form='formatted',status='unknown')
					!open (unit=2000,file='sphr_center.inp',form='formatted',status='unknown')
					write (2000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ', ' "UX" '
					do m=1,nbody
						write (2000,'(4(2x,f12.8))') xc(m,1), xc(m,2), xc(m,3), radbdy(m)
					enddo
					close(2000, status="keep")
#endif
				endif
				call alloc_phase_related
				phase_array(1)%npart = nbody
				final_vol_frac = nbody*fourthirdpi*(rad_tmp**3)/(ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
#if 0
				if (final_vol_frac.lt.matern_treshold*phiavg) then 
					call dealloc_bnd_related
					call dealloc_phase_related

					print*,'REdoING MATERN'
					goto 2002
				endif
#endif

				if (I_AM_NODE_ZERO) then
					print*,'number of bodies after scaling= ', nbody
					write (*,*)'final volume fraction after Matern: ', final_vol_frac
					write (ounit,*)'number of bodies after scaling = ', nbody

					!NOw calculate the gof 
					!call calc_gofr(nbody, xc(1:nbody,1:3), radbdy(1:nbody), mytmp, mbox,&
					!     & xperiodic,nrbins, 3, rescaling, gofr(1:nrbins),&
					!     & rho_est(1:nrbins), rad_bin(1:nrbins)) 

					gofr_avg = zero 
					do j=1,nbins
						gofr_avg(j) = gofr_avg(j) + one/dble(nsim)*sum(gofr_mis(1:nsim,j))
					enddo

#if 0
					do j=1,nbins
						if (nsim.ge.2) then 
							conf=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* sum((gofr_mis(1:nsim,j)- gofr_avg(j))**2))
						endif

						write (gofavgunit,'(10(E20.10,1x))') rad_bin(j), rad_bin(j)*lybyd, gofr_avg(j), conf
					enddo
					close(gofavgunit, status = "keep") 

					int_dist_avg = sum(int_dist(1:nsim))/nsim
					if (nsim>2) then
						call get_confin(nsim, confint)
						int_dist_var = sum( (int_dist(1:nsim)- int_dist_avg) **2) / nsim
						int_dist_sd  = sqrt(int_dist_var)

						conf = int_dist_sd / sqrt(float(nsim)) * confint
					else
						conf = zero
					endif
#endif

					!call int_dist_gofr(nbins, rad_bin(1:nbins), gofr_avg(1:nbins), int_dist2)
					!filename = "NMIS_interparticle.dat"
					!open (unit=1,file=trim(filename),status="replace",action="write")
					!write (1,"(4D15.7)") final_vol_frac, int_dist_avg, conf, int_dist2
					!close (1)

					if (GOF_AVG) then 
						if (I_AM_NODE_ZERO)write (*,*) 'GOF AVG IS true, SO stopPING THE SIMULATION'
						PARALLEL_FINISH()
						stop 
					endif
				endif
				!stop
			elseif (trim(input_type).eq."risers") then 
				doml(2) = widthbyd*dia_phys
				write (*,*)'doING RISER SIMULATION'
				if (my.eq.undefined_I) then
					my = widthbyd*dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					mz = my
					!my2 =  my/2+1
					!mx2 = my2
					if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY
				endif
				if (mx.eq.undefined_I) then
					mx = aspect_ratio*my !new + 1
					!mxf = mx
					!mx1 = mx !- 1 
					!mx2 = mx/2+1

					if (I_AM_NODE_ZERO)print*,'MX UNdeFINED: SETTING TO:', MX
				endif
				dx = doml(2)/(my)
				dy = dx
				dz = dx 

				!new doml(1) = (mx-1)*dx
				doml(1) = mx*dx
				doml(3) = mz*dz

				nphases = 1
				call alloc_phase_related

				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = volume_fraction
				char_length = dia_phys

				rad_tmp = half*(my)/lybyd
				vol_tmp =  volume_fraction
				lambda_p = (-one/(fourthirdpi*(hhat**3)*(rad_tmp**3)))
				if (one.gt.(hhat**3.0)*volume_fraction) then 
					lambda_p = lambda_p*log(one-(hhat**3.0)*volume_fraction)
				else
					write (*,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'
					write (ounit,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'
					PARALLEL_FINISH()
					stop
				endif

				if (I_AM_NODE_ZERO) then
					!new ldomtmp(1) = mx1 !-1
					ldomtmp(1) = mx
					ldomtmp(2:3) = my
					!!$ldomtmp(1) = doml(1)
					!!$ldomtmp(2:3) = doml(2:3)

					nbdymax = 2*Nint(lambda_p*ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
					write (*,*)'MAX Number of bodies in Matern : ', nbdymax
					if (.not.allocated(xc_gener))allocate(xc_gener(nbdymax,3), rad_gener(nbdymax))

20020				continue             
					call matern(3, ibordel,Ldomtmp,lambda_p, rad_tmp,hhat,mis_mat ,nbdymax,xc_gener,nbody)!0.006 
					final_vol_frac = nbody*fourthirdpi*(rad_tmp**3)/(ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
					if (final_vol_frac.lt.0.9*phiavg) then 
						print*,'REdoING MATERN'
						goto 20020
					endif
					write (*,*)'number of bodies after Matern= ', nbody
					write (*,*)'final volume fraction after Matern: ', final_vol_frac
#if 0             
					percent_buf(1:nphases) = zero
					min_part_sep = zero
					allocate(npart(nphases))
					npart(1) = nbody
					write (*,*)'nbody: ', nbody
					write (*,*)'doml(:) ', doml(:)
					read (*,*)

					do IDIM = 1, 3
						XC_GENER(1:NBODY,IDIM) = XC_GENER(1:NBODY, IDIM)/doml(IDIM)
					enddo

					rad_gener(1:NBODY)=rad_tmp/doml(2)

					call scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
						&,percent_buf(1:nphases),xc_gener(1:nbody,1:3),rad_gener(1:nbody),&
						& min_part_sep, toscale=.true.)
					write (*,*),'nbody: ', nbody
#endif                  
					call alloc_bnd_related
					do n = 1, ndim
						xc(1:nbody,n) = xc_gener(1:nbody,n) + one
					enddo
					radbdy(1:nbody) = rad_tmp
				endif

				do iphs = 1, nphases
					phase_array(iphs)%npart = nbody
				enddo
          
				BROADCAST_INT(nbody,1,node_zero,comm_cart_2d)

				if (.not.I_AM_NODE_ZERO)call alloc_bnd_related

				do idim = 1, ndim
					BROADCAST_DOUBLE(xc(1,idim),nbody,node_zero,comm_cart_2d)
				enddo
				BROADCAST_DOUBLE(radbdy(1),nbody,node_zero,comm_cart_2d)
			elseif (input_type.eq."random") then 
				call generate_psd_config
				do idim = 1, ndim
					BROADCAST_DOUBLE(xc(1,idim),nbody,node_zero,comm_cart_2d)
				enddo
				BROADCAST_DOUBLE(radbdy(1),nbody,node_zero,comm_cart_2d)
			elseif (input_type.eq."simple") then
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'S')
					write (*,*) 'IN SIMPLE CUBIC'
				endif
				nbody = 1
				nphases = 1
				call alloc_phase_related
				phase_array(1)%npart = nbody
				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = phiavg
				char_length = dia_phys
				call alloc_bnd_related

				lybyd = (pi/(6.d0*phiavg))**(one/three)
				doml(2) = lybyd*dia_phys
				if (my.eq.undefined_I) then
					my = lybyd * dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					!if (my/2.ne.0)my = my + 1
					mz = my
					mx = my !new +1
					!mxf = mx 
					!my2 =  my/2+1
					!mx1 = mx !new- 1 
					!mx2 = mx/2+1

					if (I_AM_NODE_ZERO) print*,'MY UNdeFINED: SETTING TO:', MY
				endif

				xc(1:nbody,1) = dble(mx-1)/two + one !+one
				xc(1:nbody,2) = dble(my-1)/two + one !+one
				xc(1:nbody,3) = dble(mz-1)/two + one !+one

				!xc(1:nbody,1) = mx/two + one
				!xc(1:nbody,2) = my/two + one
				!xc(1:nbody,3) = mz/two + one

				!xc(1,1:3) = my/2.d0+one
				xperiodic = .true.
       
				dx = doml(2)/(my)
				dy = dx
				dz = dx 
				!new doml(1) = (mx-1)*dx
				doml(1) = mx*dx
				doml(3) = mz*dz

				radbdy(1:nbody) = half*(my)/lybyd
				!xc(1:nbody, 1) =  radbdy(1:nbody)!mx1/two+one
				if (I_AM_NODE_ZERO) then
					print*,'phiavg: ', phiavg
					!print*,'phiavg in SIMPLE ACTUAL: ', nbody*
					print*,'LYBYD AND RADBDY',lybyd, radbdy(1)
					call screen_separator(80,'S')
				endif
			elseif (input_type.eq."fcc") then
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'FCC')
					write (*,*) 'IN FACE CENTERED CUBIC'
				endif
				nbody = 4
				nphases = 1
				call alloc_phase_related
				phase_array(1)%npart = nbody
				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = phiavg
				char_length = dia_phys
				call alloc_bnd_related
       
				lybyd= (4.d0*pi/(6.d0*phiavg))**(one/three)
				if (my.eq.undefined_I) then
					my = lybyd * dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					!if (my/2.ne.0)my = my + 1
					mz = my
					mx = my !new +1
					!mxf = mx 
					!my2 =  my/2+1
					!mx1 = mx !new - 1 
					!mx2 = mx/2+1

					if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY
				endif

				xc(1,1:3) = 1.d0
				xc(2,1) = my/2.d0+one
				xc(2,2) = my/2.d0+one
				xc(2,3) = 1.d0
				xc(3,1) = my/2.d0+one
				xc(3,2) = 1.d0
				xc(3,3) = my/2.d0+one
				xc(4,1) = 1.d0
				xc(4,2) = my/2.d0+one
				xc(4,3) = my/2.d0+one

				xperiodic = .true.

				doml(2) = lybyd*dia_phys
				dx = doml(2)/(my)
				dy = dx
				dz = dx 
				!new doml(1) = (mx-1)*dx
				doml(1) = mx*dx
				doml(3) = mz*dz
       
				radbdy(1:nbody) = half*(my)/lybyd
				if (I_AM_NODE_ZERO) then
					print*,'phiavg: ', phiavg
					!print*,'phiavg in SIMPLE ACTUAL: ', nbody*
					print*,'LYBYD AND RADBDY',lybyd, radbdy(1)
					call screen_separator(80,'FCC')
				endif
			elseif (input_type.eq."default") then
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'D')
					write (*,*) 'IN DEFAULT FRESH START'

					tmp_unit =  getnewunit(minunitno,maxunitno)
					open (tmp_unit,file=trim(run_name)//'_sphr_center.inp',status='old', Action='read')

					phiavg = zero
					nbody = 0
					!read (tmp_unit,*)
125				        continue
					read (tmp_unit,*,end=225) dummyx, dummyy, dummyz, dummyr ! , dummyi !M ,dummyx,dummyy,dummyz,dummyx
					!if (dummyy<=lybyd/8.or.dummyy>=lybyd*7/8) then
					if (dummyx<=lybyd.and.dummyy<=lybyd.and.dummyz<=lybyd) nbody = nbody + 1
					!endif
					goto 125
225				        continue
					close(tmp_unit, status='keep')
					write (*,*)'Number of bodies in the deFAULT CASE are : ', nbody
				endif
				BROADCAST_INT(nbody,1,node_zero,comm_cart_2d)
				if (nbody.NE.0) then
					call alloc_bnd_related
					if (I_AM_NODE_ZERO) then
						open (tmp_unit,file=trim(run_name)//'_sphr_center.inp',status='old', Action='read')
						write (*,*) 'Reading immersed surface data'
						!read (tmp_unit,*)
						i=0
						do
							read (tmp_unit,*) dummyx, dummyy, dummyz, dummyr !, dummyi

							!if (dummyy<=lybyd/8.or.dummyy>=lybyd*7/8) then
							if (dummyx<=lybyd.and.dummyy<=lybyd.and.dummyz<=lybyd) then 
								i=i+1
								xc(i,1) = dummyx*dbydx
								xc(i,2) = dummyy*dbydx
								xc(i,3) = dummyz*dbydx
								!color(i) = dble(dummyi)

								velbdy(i,1:ndim) = zero
								!M read (tmp_unit,*) xc(i,1),xc(i,2),xc(i,3), radbdy(i), elbdy(i,1), velbdy(i,2), velbdy(i,3)
								!M ol_frac1 = phiavg + (pi*(dia_phys**3.d0))/(6.d0)
							endif
							!endif
							if (i==nbody) exit
						enddo
						if (i/=nbody) then
							write (*,*) "ERROR IN READING INPUT DATA..."
							stop
						else
							write (*,*)'Number of bodies in the deFAULT CASE are : ', nbody
						endif


						close(tmp_unit,STATUS='keep')
						!dummyx = MINVAL(xc(1:nbody,1))
						!do m = 1, nbody
						!	xc(m, 1) = xc(m, 1) - dummyx + 1
						!enddo
					endif
					nphases = 1
					call alloc_phase_related
					phase_array(1)%npart = nbody
					phase_array(1)%dia = dia_phys
					phase_array(1)%volfrac = phiavg
					char_length = dia_phys
				endif
		    
				do idim = 1, ndim
					BROADCAST_DOUBLE(xc(1,idim),nbody,node_zero,comm_cart_2d)
					BROADCAST_DOUBLE(velbdy(1,idim),nbody,node_zero,comm_cart_2d)
				enddo
#if 1
				allocate(npart(nphases))
				do iphs = 1, nphases
					npart(iphs) = nbody
					percent_buf(iphs) = zero
				enddo
#endif             
				if (I_AM_NODE_ZERO) print*,'PART VOLUME = ', phiavg
				doml(2) = lybyd*dia_phys

				if (my.eq.undefined_I) then
					my = lybyd * dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					!if (my/2.ne.0)my = my + 1
					mz = my
					mx = my !new +1
					!mxf = mx 
					!my2 =  mx/2+1
					!mx1 = mx !new - 1 
					!mx2 = mx1/2+1
					if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY
				endif

				!Convert the centers to lie between 0 and 1
				do m = 1, nbody
					! xc(m,1:ndim) = xc(m,1:ndim) * real(my,prcn)/(lybyd*dbydx)
					radbdy(m) = half*real(my,prcn)/lybyd !M half*one/lybyd
				enddo

				dx = doml(2)/(my)
				dy = dx
				dz = dx 
				!new doml(1) = (mx-1)*dx
				doml(1) = mx*dx
				doml(3) = mz*dz
#if 0
				call scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
					&,percent_buf(1:nphases),xc(1:nbody,1:3),radbdy(1:nbody),&
					& min_part_sep, toscale=.true.) 

				do iphs = 1, nphases
					phase_array(iphs)%npart = npart(iphs)
				enddo
#endif             !stop
       
				if (I_AM_NODE_ZERO) then
					write (*,*) 'OUT OF deFAULT SPHR CONFIG'
					call screen_separator(80,'D')
				endif
			elseif (input_type.eq."agglomerate") then
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'A')
					write (*,*) 'INPUT TYPE AGGLOMERATE FRESH START'

					tmp_unit =  getnewunit(minunitno,maxunitno)
					open (tmp_unit,file=trim(run_name)//'_sphr_center.inp',status='old', Action='read')

					phiavg = zero
					nbody = 0
					!read (tmp_unit,*)
1250				continue
					read (tmp_unit,*,end=2250) dummyx, dummyy, dummyz
					nbody = nbody + 1
					goto 1250
2250				continue
					close(tmp_unit, status='keep')
					write (*,*)'Number of bodies in the deFAULT CASE are : ', nbody
				endif
				BROADCAST_INT(nbody,1,node_zero,comm_cart_2d)

				doml(2) = lybyd*dia_phys
				if (my.eq.undefined_I) then
					my = lybyd * dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					!if (my/2.ne.0)my = my + 1
					mz = my
					mx = my !new +1
					!mxf = mx 
					!my2 =  my/2+1
					!mx1 = mx !new - 1 
					!mx2 = mx/2+1

					if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY
				endif
				dx = doml(2)/(my)
				dy = dx
				dz = dx 
				!new doml(1) = (mx-1)*dx
				doml(1) = mx*dx
				doml(3) = mz*dz

				if (nbody.NE.0) then
					call alloc_bnd_related
					if (I_AM_NODE_ZERO) then
						open (tmp_unit,file=trim(run_name)//'_sphr_center.inp',status='old', Action='read')
						write (*,*) 'Reading immersed surface data'
						!read (tmp_unit,*)
						do i = 1, nbody
							read (tmp_unit,*) xc(i,1),xc(i,2),xc(i,3)
							velbdy(i,1:ndim) = zero
						enddo

						close(tmp_unit,STATUS='keep')

						do m = 1, nbody
							do idim = 1, ndim
								xc(m, idim) = xc(m, idim)/doml(idim) + half
							enddo
						enddo
					endif
					nphases = 1
					call alloc_phase_related
					phase_array(1)%npart = nbody
					phase_array(1)%dia = dia_phys
					phase_array(1)%volfrac = phiavg
					char_length = dia_phys
				endif
          
#if PARALLEL
				do idim = 1, ndim
					BROADCAST_DOUBLE(xc(1,idim),nbody,node_zero,comm_cart_2d)
					BROADCAST_DOUBLE(velbdy(1,idim),nbody,node_zero,comm_cart_2d)
				enddo
#endif

				allocate(npart(nphases))
				do iphs = 1, nphases
					npart(iphs) = nbody
					percent_buf(iphs) = zero
				enddo
	       
				do m = 1, nbody
					radbdy(m) = half*(one/lybyd)
				enddo
	       
				call scale_to_grid_units(nbody,npart(1:nphases),nphases,my,my,xperiodic(1)&
					&,percent_buf(1:nphases),xc(1:nbody,1:3),radbdy(1:nbody),&
					& min_part_sep) 

				do iphs = 1, nphases
					phase_array(iphs)%npart = npart(iphs)
				enddo
	       
				if (I_AM_NODE_ZERO) then
					write (*,*) 'OUT OF AGGLOMERATE FRESH START'
					call screen_separator(80,'A')
				endif
			elseif (trim(input_type).eq.'lubtest') then
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'L')
					write (*,*) 'IN LUBRICATION THEORY TEST --> FRESH START'
					!!$ tmp_unit =  getnewunit(minunitno,maxunitno)
					!!$ open (tmp_unit,file=trim(run_name)//'_lubtest.inp',status='old', Action='read')
					!!$ read (tmp_unit,*)nphases
				endif
				nphases = 2
				call alloc_phase_related
				allocate(npart(nphases))
				do iphs = 1, nphases
					phase_array(iphs)%npart = 1
					npart(iphs) = 1
					percent_buf(iphs) = zero
				enddo
				phase_array(1)%dia = dia_phys
				phase_array(2)%dia = phase_array(1)%dia*dia_ratio
		       
				char_length = two*(phase_array(1)%dia*phase_array(2)%dia)/(phase_array(1)%dia+phase_array(2)%dia)

				doml(2) = lybyd*phase_array(2)%dia

				doml(1) = doml(2)
				doml(3) = doml(2)

				nbody = 2
				call alloc_bnd_related

				do m = 1, nbody
					radbdy(m) = phase_array(m)%dia/two
					xc(m,2) = doml(2)/two
					xc(m,3) = doml(2)/two
				enddo
				part_sep = hbyd*phase_array(2)%dia

				xc(1,1) = (doml(2)/two - part_sep/two - radbdy(1))
				xc(2,1) = (doml(2)/two + part_sep/two + radbdy(2))

				do m = 1, nbody
					radbdy(m) = radbdy(m)/doml(2)
					do idim = 1, ndim
						xc(m,idim) = xc(m,idim)/doml(idim)
					enddo
				enddo
		       
				if (I_AM_NODE_ZERO) then
					write (*,'(A25,2x,g17.8)')'LENGTH OF THE BOX IS',doml(2)
					write (*,'(A25,2(2x,g17.8))')'RADII ',radbdy(1), radbdy(2)
					write (*,'(A25,2(2x,g17.8))')'PHYSICAL LOCATIONS ',xc(1,1), xc(2,1)
				endif
          
				if (my.eq.undefined_I) then
					my = doml(2)/phase_array(1)%dia * dbydx
					if (MOD(my,2).ne.0) then
						my = my+1
						print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
					endif
					!if (my/2.ne.0)my = my + 1
					mz = my
					mx = my !new +1
					!mxf = mx 
					!my2 = my/2+1
					!mx1 = mx !enw - 1 
					!mx2 = mx/2+1
					if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY
				endif
				dx = doml(2)/real(my,prcn)
				dy = dx
				dz = dx

				!if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY

				min_part_sep = zero
				hbydx = part_sep/dx
				write (*,'(A25,2(2x,g17.8))')'SEPARATION IN GRID UNITS : ',hbydx
				call scale_to_grid_units(nbody,npart(1:nphases),nphases,my,my,xperiodic(1)&
					&,percent_buf(1:nphases),xc(1:nbody,1:3),radbdy(1:nbody),&
					& min_part_sep) 
		       
				if (I_AM_NODE_ZERO) then
					print*,'RADII IN GRID UNITS: ',radbdy(1), radbdy(2)
					call screen_separator(80,'L')
				endif
			elseif (trim(input_type).eq.'single-phase') then
				doml(1:ndim) = two*pi
!				call screen_separator(80,'T')
!				write (*,*) 'IN SINGLE PHASE TURBULENCE RUN --> FRESH START'
!				if (my.eq.undefined_I) then
!					write (*,'(A)')'MY IS UNdeFINED FOR SINGLE-PHASE RUN. deFINE IT In floparam.in AND TRY AGAIN.'
!					PARALLEL_FINISH()
!					stop
!				endif
				mz = my
				mx = my !new +1

				dx = doml(1)/real(mx,prcn)
				dy = doml(2)/real(my,prcn)
				dz = doml(3)/real(mz,prcn)

				nbody = 0
				nphases = 0
				char_length = doml(2)
			endif !ihardsphere 
		else !irestart
			if (I_AM_NODE_ZERO) then
				open (unit=1050,file=trim(run_name)//'_RESTART',form="formatted",status="unknown", action="read")    
				read (1050, *) count_restart, nproc, nprocy, nprocz
				close(1050, status='keep')

				if (imove==1 .and. .not.moving_from_fixed_bed) then
					write (filename1,'(I1)') count_restart
					open (1050,file=trim(run_name)//'_sphere_config_'//trim(filename1)//'.rst',form=rstsave, status="unknown")
				else
					open (1050, file = trim(run_name)//'_sphere_config.rst', form="unformatted", status="old", action="read")
				endif

				nbody = 0
				nphases = 0
				read (1050) nbody
				if (nbody>0) then
					read (1050) nphases
					call alloc_phase_related
					do iphs = 1, nphases
						read (1050) phase_array(iphs)%npart
					enddo
					do iphs = 1, nphases
						read (1050) phase_array(iphs)%dia
					enddo
					read (1050) doml(1:ndim)
					read (1050) mx,my,mz
					read (1050) dx,dy,dz
					dy = doml(2)/real(my,prcn)
					dx = dy
					dz = dy 
					read (1050) MOVE_PARTICLES
					read (1050) char_length
				elseif (trim(input_type).eq."lubtest") then
					nphases = 2
					call alloc_phase_related
				elseif (trim(input_type).eq."single-phase") then
					nphases = 0

					read (1050) doml(1:ndim)
					read (1050) mx,my,mz
					read (1050) dx,dy,dz
					dy = doml(2)/real(my,prcn)
					dx = dy
					dz = dy
					read (1050) char_length
				else
					nphases = 1
					call alloc_phase_related
				endif
				write (*,*) 'FROM RESTART: NBODY = ', nbody
				write (*,"(1a,3d15.7)") 'FROM RESTART: DOML     = ', doml(1:ndim)
				write (*,"(1a,3d15.7)") 'FROM RESTART: DX,DY,DZ = ', dx, dy, dz
				write (*,"(1a,3i)")     'FROM RESTART: MX,MY,MZ = ', mx, my, mz
			endif

			BROADCAST_INT(count_restart,1,node_zero,comm_cart_2d)
			BROADCAST_INT(nbody,1,node_zero,comm_cart_2d)
			BROADCAST_INT(nphases,1,node_zero,comm_cart_2d)

			if (.not.I_AM_NODE_ZERO) then
				if (nphases.gt.0) call alloc_phase_related
			endif

			if (nphases>0) then
				do iphs=1, nphases
					BROADCAST_INT(phase_array(iphs)%npart,1,node_zero,comm_cart_2d)
					BROADCAST_DOUBLE(phase_array(iphs)%dia,1,node_zero,comm_cart_2d)
				enddo
			endif
			BROADCAST_DOUBLE(doml, ndim, node_zero, comm_cart_2d)
			BROADCAST_INT(mx, 1, node_zero, comm_cart_2d)
			BROADCAST_INT(my, 1, node_zero, comm_cart_2d)
			BROADCAST_INT(mz, 1, node_zero, comm_cart_2d)
			BROADCAST_DOUBLE(dx, 1, node_zero, comm_cart_2d)
			BROADCAST_DOUBLE(dy, 1, node_zero, comm_cart_2d)
			BROADCAST_DOUBLE(dz, 1, node_zero, comm_cart_2d)
			BROADCAST_LOGICAL(MOVE_PARTICLES, 1, node_zero, comm_cart_2d)
			BROADCAST_DOUBLE(char_length, 1, node_zero, comm_cart_2d)

			if (nbody>0) then
				call alloc_bnd_related
				if (I_AM_NODE_ZERO) then
					read (1050) xc(1:nbody,1:3)
					read (1050) radbdy(1:nbody)
					read (1050) velbdy(1:nbody,1:3)

					read (1050) frame_vel(1:ndim)
					read (1050) frame_pos(1:ndim)
					close(1050, status = "keep")
					!------------------------------------------------------------------------
				endif
			endif

			if (nbody>0) then
				do idim = 1, ndim
					BROADCAST_DOUBLE(xc(1,idim),nbody,node_zero,comm_cart_2d)
					BROADCAST_DOUBLE(velbdy(1,idim),nbody,node_zero,comm_cart_2d)
				enddo
				BROADCAST_DOUBLE(radbdy(1),nbody,node_zero,comm_cart_2d)

				BROADCAST_DOUBLE(frame_vel,3,node_zero,comm_cart_2d)
				BROADCAST_DOUBLE(frame_pos,3,node_zero,comm_cart_2d)
			endif

			if (input_type.eq."random".or.trim(input_type).eq.'risers') then
				do iphs = 1, nphases
					BROADCAST_INT(phase_array(iphs)%npart, 1, node_zero, comm_cart_2d)
					BROADCAST_DOUBLE(phase_array(iphs)%dia, 1, node_zero, comm_cart_2d)
				enddo
			elseif (input_type.eq."simple") then
				xperiodic = .true.
				phase_array(1)%npart = nbody
				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = phiavg
				lybyd= (pi/(6.d0*phiavg))**(one/three)
				if (I_AM_NODE_ZERO)print*,'LYBYD IN SIMPLE',lybyd
				!doml(2) = lybyd*dia_phys
			elseif (input_type.eq."fcc") then
				xperiodic = .true.
				lybyd= (4.d0*pi/(6.d0*phiavg))**(one/three)
				phase_array(1)%npart = nbody
				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = phiavg
				if (I_AM_NODE_ZERO) print*,'LYBYD IN FCC',lybyd
				!doml(2) = lybyd*dia_phys
			elseif (input_type.eq."default") then 
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'D')
					write (*,*) 'IN deFAULT SPHR CONFIG RESTART'
				endif
				phiavg = zero 
				do i = 1, nbody
					phiavg = phiavg + (pi*(dia_phys**3.d0))/(6.d0)
#if 0
					if (I_AM_NODE_ZERO) then
						write (*,*) 'NBODY = ', NBODY
						write (*,*) 'XC = ', XC(1,:)
						write (*,*) 'RADBDY and VELBDY =', RADBDY(I), VELBDY(I,:)
					endif
#endif
				enddo
				phase_array(1)%npart = nbody
				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = phiavg

				!doml(2) = lybyd*dia_phys

				if (I_AM_NODE_ZERO) print*,'PART VOL =  ', phiavg

				if (I_AM_NODE_ZERO) then
					write (*,*) 'OUT OF deFAULT SPHR CONFIG RESTART'
					call screen_separator(80,'D')
				endif
			elseif (input_type.eq."agglomerate") then 
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'A')
					write (*,*) 'IN AGGLOMERATE RESTART'
				endif
				phiavg = zero 
				do i = 1, nbody
					phiavg = phiavg + (pi*(dia_phys**3.d0))/(6.d0)
#if 0
					if (I_AM_NODE_ZERO) then
						write (*,*) 'NBODY = ', NBODY
						write (*,*) 'XC = ', XC(1,:)
						write (*,*) 'RADBDY and VELBDY =', RADBDY(I), VELBDY(I,:)
					endif
#endif
				enddo
				phase_array(1)%npart = nbody
				phase_array(1)%dia = dia_phys
				phase_array(1)%volfrac = phiavg
				!doml(2) = lybyd*dia_phys
				print*,'PART VOL =  ', phiavg

				if (I_AM_NODE_ZERO) then
					write (*,*) 'OUT OF AGGLOMERATE RESTART'
					call screen_separator(80,'A')
				endif
			elseif (input_type.eq."lubtest") then 
				if (I_AM_NODE_ZERO) then
					call screen_separator(80,'L')
					write (*,*) 'IN LUBRICATION TEST RESTART'
				endif

				do i = 1, nbody
					if (I_AM_NODE_ZERO) then
						write (*,*) 'NBODY = ', NBODY
						write (*,*) 'XC = ', XC(1,:)
						write (*,*) 'RADBDY and VELBDY =', RADBDY(I), VELBDY(I,:)
					endif
				enddo
				do iphs = 1, nphases
					phase_array(iphs)%npart = 1
				enddo
				phase_array(1)%dia = dia_phys
				phase_array(2)%dia = phase_array(1)%dia*dia_ratio
				!doml(2) = lybyd*phase_array(2)%dia

				if (I_AM_NODE_ZERO) then
					write (*,*) 'OUT OF LUBRICATION TEST RESTART'
					call screen_separator(80,'L')
				endif
			elseif (trim(input_type).eq.'single-phase') then
				doml(1:ndim) = two*pi
				char_length = doml(2)
				call screen_separator(80,'T')
				write (*,*) 'IN SINGLE PHASE TURBULENCE RUN --> RESTART'
				if (my.eq.undefined_I) then
					write (*,'(A)')'MY IS UNdeFINED FOR SINGLE-PHASE RUN. deFINE IT In floparam.in AND TRY AGAIN.'
					PARALLEL_FINISH()
					stop
				endif
				mz = my
				mx = my !new +1
				!mxf = mx 
				!my2 =  my/2+1
				!mx1 = mx !new - 1 
				!mx2 = mx/2+1
			endif
       
			if (my.eq.undefined_I) then
				!if ((trim(input_type).eq.'random').and.(trim(psd_type).eq.'bidisp')) then
				if ((trim(input_type).eq.'random').and.(trim(psd_type).eq.'csd')) then !Need to be done correctly for CSD
					my = doml(2)/phase_array(1)%dia * dbydx
				elseif (trim(input_type).eq.'risers') then
					my = widthbyd*dbydx
				elseif (trim(input_type).eq.'lubtest') then
					my = doml(2)/phase_array(1)%dia * dbydx
				else
					my = lybyd * dbydx
				endif
				if (MOD(my,2).ne.0) then
					my = my+1
					print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
				endif
				!if (my/2.ne.0)my = my + 1
				mz = my
				!my2 =  my/2+1
				if (I_AM_NODE_ZERO)print*,'MY UNdeFINED: SETTING TO:', MY
			endif
			if (mx.eq.undefined_I) then
				if (trim(input_type).eq.'risers') then
					mx = aspect_ratio*my !new + 1
				else
					mx = my !new +1
				endif
				!mxf = mx 
				!mx1 = mx !new - 1
				!mx2 = mx/2
			endif
       
			dy = doml(2)/my
			dy = dx
			dz = dx
			!new doml(1) = (mx-1)*dx
			doml(1) = mx*dx
			doml(3) = mz*dz

#if 0
			if (nbody>0) then
				!nbins = 200
				if (.not.allocated(gofr_avg)) then
					nsim = 1
					allocate(gofr_avg(nbins), gofr_mis(nsim, nbins))
					allocate(rad_bin(nbins))
					allocate(int_dist(nsim))

					gofr_avg = zero
					gofr_mis = zero
					rad_bin  = zero
					int_dist = zero
				endif

				write (*,*) 'GENERATING GOFR' 
				if (.not.allocated(contact)) allocate(contact(nbody,nbody))
				call calculate_gofr_homog(nbody,xc(1:nbody,1:3), contact, global_n(:), nbins, .true., gofr(1:nbins), rad_bin(1:nbins), max_overlap)

				do j = 1, nsim
					gofr_mis(j,:) = gofr(:)
				enddo

				gofr_avg = zero 
				do j=1,nbins
					gofr_avg(j) = gofr_avg(j) + one/dble(nsim)*sum(gofr_mis(1:nsim,j))
				enddo
			  
				gofavgunit = getnewunit(minunitno,maxunitno)
				open (gofavgunit,file= trim(run_name)//'_gof_avg.dat', form = 'formatted')


				do j=1,nbins
					write (gofavgunit,'(3d15.7)') rad_bin(j), rad_bin(j)*lybyd, gofr_avg(j)
				enddo
				close(gofavgunit, status = "keep") 
			endif
#endif       
		endif !irestart
    
#if 0
		if (I_AM_NODE_ZERO.and.nbody.ne.0) then
			open (unit=2000,file=trim(run_name)//'_sphr_center_out.dat',form='formatted',status='unknown')
			write (2000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" ', ' "IBDY" '
			!!$ write (2000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ', ' "UX" ', ' "UY" ', ' "UZ" ', ' "PHI" '

			do m=1,nbody
				!write (2000,'(4(2x,g17.8),2x,i3)') xc(m,1), xc(m,2), xc(m,3), radbdy(m),m
				write (2000,'(6(2x,g17.8))') xc(m,1), xc(m,2), xc(m,3), radbdy(m)
			enddo
			close(2000, status="keep")
		endif
#endif

		if (I_AM_NODE_ZERO) write (*,*)'DONE WITH RESTART'

	end subroutine generate_configuration

	subroutine generate_particle_velocities
		use global_data
		use mypost_process, only : uf_cage, fluid_vel_cage

		implicit none

		real(prcn) :: lambda_p, gran_energy_out,&
			& diff_vel_spec1(ndim), diff_vel_spec2(ndim), temp_spec1,&
			& temp_spec2, mean_gt, umean_test1(ndim), umean_test2(ndim),&
			& mean_vel(ndim),gran_energy_out1,gran_energy_out2,&
			& var_test(ndim), hist(20), fmin, fmax, umf0(ndim),&
			& rsf(ndim,ndim), rsf1(ndim,ndim), rsf2(ndim,ndim), x1, x2,&
			& vfrac_mean,variance1,variance2
		real(prcn) :: gran_temp_spec(nphases),diff_vel_spec(nphases,ndim)&
			&, constant1(ndim), constant2, test_mix_var,&
			& test_mixmean_vmag, test_mean_spec_vel(nphases,ndim),&
			& test_mixmean_vel(ndim), test_spec_var(nphases)

		real(prcn), dimension(:,:), allocatable ::  vtemp,vtemp1,vtemp2,veltemp
		real(prcn), dimension(:), allocatable ::  ftemp, wt 
		integer ::  unit1, ncount(ndim), count_temp(ndim),idim, i, j, iphs, pstart, pend, m

		velbdy(1:nbody,:) = zero

		!if (nphases.gt.1)mean_vel_to_particles = .true.
		!Assign Mean velocity to the particles

		if (init_parts_with_fluid_vel) then
			call fluid_vel_cage

			do m=1, nbody
				velbdy(m,:) = uf_cage(m,:)
			enddo
			deallocate(uf_cage)
			goto 10
		endif

		if (mean_vel_to_particles) then
			if (trim(psd_type).eq."discrete".or.trim(psd_type).eq."bidisp") then
				if (equal_momentum) then
					constant2 = zero
					do iphs = 1, nphases
						constant2 = constant2 + phase_array(iphs)%npart/phase_array(iphs)%volfracg
					enddo
					constant2 = constant2*maxvolfrac/((one-maxvolfrac)*real(nbody,prcn)) + real(nphases,prcn)/maxvolfrac
					do idim = 1, ndim
						!!$constant1(idim) = maxvolfrac*uchar(idim)/((one-maxvolfrac)*real(nphases,prcn))
						constant1(idim) = uchar(idim)/((one-maxvolfrac)*(constant2))                   
					enddo
					do iphs = 1, nphases
						do idim = 1, ndim
							phase_array(iphs)%mean_spec_vel(idim) = constant1(idim)/phase_array(iphs)%volfracg
						enddo
					enddo
				else
					do iphs = 1, nphases
						do idim = 1, ndim
							phase_array(iphs)%mean_spec_vel(idim) = -uchar(idim)/(one-maxvolfrac) ! Equal velocities
							!!$phase_array(iphs)%mean_spec_vel(idim) = uchar(idim) ! Equal velocities
						enddo
					enddo
				endif
				pstart = 1
				do iphs = 1, nphases
					do idim=1,ndim
						diff_vel_spec(iphs,idim) = phase_array(iphs)%mean_spec_vel(idim)- uchar(idim)
					enddo
				enddo
         
				do m = 1, nbody
					iphs = part_array(m)%iphs
					velbdy(m,1:ndim) = phase_array(iphs)%mean_spec_vel(1:ndim)
				enddo
			elseif (psd_type.eq."csd") then
				do idim = 1, ndim  
					velbdy(1:nbody,idim) = uchar(idim)!/(one-maxvolfrac) ! Assigning equal velocities to all the spheres for the time being. Need to decide on Re_m(r).
				enddo
			endif
		else
			do iphs = 1, nphases
				do idim = 1, ndim
					phase_array(iphs)%mean_spec_vel(idim) = zero
				enddo
			enddo
			velbdy(1:nbody,1:ndim) = zero
		endif

		! Generate velocity distribution to the particles
		if (ReT.gt.SMALL_NUMBER) then
			gran_temp = (ReT*vis/dia_phys)**2.d0
			if (I_AM_NODE_ZERO) then
				if (psd_type.eq."csd") then
					do j=1,ndim
						umf0(j)=0.0
						do i=1,ndim
							if (i.eq.j) then
								rsf(i,j)= gran_temp
							else
								rsf(i,j)=0.0
							endif
						enddo
					enddo
					allocate(vtemp(nbody,ndim))
					call jn_dist(vtemp,SIZE(vtemp,1),ndim,umf0,rsf) ! Maxwellian
					! distribution for the particle velocity fluctuations
					do m = 1, nbody
						velbdy(m,1:ndim) = velbdy(m,1:ndim) + vtemp(m,1:ndim)
					enddo
					deallocate(vtemp)
				elseif (trim(psd_type).eq."discrete".or.trim(psd_type).eq."bidisp") then
					!!$ constant2 = three*SUM(volfrac(1:nphases))*gran_temp
					constant2 = mean_volfrac*gran_temp
					do iphs = 1, nphases
						constant2 = constant2 !-  doT_PRODUCT(diff_vel_spec(iphs,1:ndim),diff_vel_spec(iphs,1:ndim))*volfrac(iphs)
					enddo
					constant2 = constant2/real(nphases,prcn)
					pstart = 1
					do iphs = 1, nphases
						gran_temp_spec(iphs) = constant2/(phase_array(iphs)%volfrac)
						!!$ gran_temp_spec(iphs) = constant2/(three*volfrac(iphs))
						do j=1,ndim
							umf0(j) = phase_array(iphs)%mean_spec_vel(j)
							do i=1,ndim
								if (i.eq.j) then
									rsf(i,j)= gran_temp_spec(iphs)
								else
									rsf(i,j)=0.0
								endif
							enddo
						enddo
						allocate(vtemp(phase_array(iphs)%npart,ndim))
						call jn_dist(vtemp,SIZE(vtemp,1),ndim,umf0,rsf) ! Maxwellian
						pend = pstart + phase_array(iphs)%npart - 1
						i = 1
						do m = pstart, pend
							velbdy(m,1:ndim) = vtemp(i,1:ndim)
							i = i+1
						enddo
						pstart = pend + 1
						deallocate(vtemp)
					enddo
				endif
			endif
		endif

10		continue
		BROADCAST_DOUBLE(velbdy, nbody*ndim, node_zero, comm_cart_2d)
	end subroutine generate_particle_velocities

  subroutine INPUT_CHECK
    implicit none
    real(prcn) :: gran_temp_spec(nphases),diff_vel_spec(nphases,ndim)&
         &, constant1(ndim), constant2, test_mix_var,&
         & test_mixmean_vmag, test_mean_spec_vel(nphases,ndim),&
         & test_mixmean_vel(ndim), test_spec_var(nphases),&
         & volfrac(nphases), vfracmean, vbox,&
         & dia(nphases), diag(nphases)
    integer :: pstart, pend, iphs, m, idim, npart(nphases)
    
    vbox = doml(1)*doml(2)*doml(3)
    volfrac = zero
    npart = 0
    do m = 1, nbody
       iphs = part_array(m)%iphs
       npart(iphs) = npart(iphs)+1
       dia(iphs) = two*radbdy(m)*dx
       diag(iphs) = two*radbdy(m)
       volfrac(iphs) = volfrac(iphs) + pi*(dia(iphs))**3.d0/6.d0
    enddo

    do iphs = 1, nphases
       volfrac(iphs) = volfrac(iphs)/vbox
    enddo

    vfracmean = SUM(volfrac(1:nphases))
    if (I_AM_NODE_ZERO) then
       call screen_separator(80,'I')
       if (.not.trim(psd_type).eq.'csd') then
          write (*,'((A25,2x,I6))')'Number of PHASES = ', nphases
          do iphs = 1, nphases
             write (*,'((A25,2x,I6, A5,2(2x,I8)))')'No. of Particles in phase ', iphs,&
                  & ' = ', npart(iphs),(phase_array(iphs)%npart)
          enddo
          
          do iphs = 1, nphases
             write (*,'((A25,2x,I6, 2x, A5,2(2x,g17.8)))')'Diameter of phase ', iphs, ' = ', dia(iphs), diag(iphs)
          enddo

          do iphs = 1, nphases
             write (*,'((A25,2x,I6, 2x, A5,2(2x,g17.8)))')'Volume fraction &
                  &of phase ', iphs, ' = ', volfrac(iphs)&
                  &,(phase_array(iphs)%volfracg)
          enddo
       endif
       write (*,'((A25,2(2x,g17.8)))')'Mean Volume fraction = ', vfracmean, maxvolfrac
    endif
    
    if (trim(psd_type).eq."discrete".or.trim(psd_type).eq."bidisp") then
       test_mean_spec_vel(1:nphases,1:ndim) = zero
       test_mixmean_vel = zero
       pstart = 1
       do iphs = 1, nphases
          pend = pstart + phase_array(iphs)%npart- 1
          do m = pstart, pend
             test_mean_spec_vel(iphs,1:ndim) = test_mean_spec_vel(iphs,1:ndim) + velbdy(m,1:ndim)
          enddo
          pstart = pend + 1
          test_mean_spec_vel(iphs,1:ndim) = test_mean_spec_vel(iphs,1:ndim)/real(phase_array(iphs)%npart, prcn)
          test_mixmean_vel(1:ndim) =  test_mixmean_vel(1:ndim) + test_mean_spec_vel(iphs,1:ndim)*phase_array(iphs)%volfracg/maxvolfrac
       enddo
       do idim = 1, ndim
          test_mixmean_vel(idim) = (one-maxvolfrac)*(test_mixmean_vel(idim)-&
               & ufmean_des(idim))
       enddo

       test_mixmean_vmag = DSQRT(doT_PRODUCT(test_mixmean_vel(1:ndim)&
            &,test_mixmean_vel(1:ndim)))

!!$       print*,'test_mean_mag = ', test_mixmean_vmag
!!$       print*,'test_mean_mag = ', phase_array(1)%volfracg,&
!!$            & phase_array(2)%volfracg, phase_array(1)%volfracg&
!!$            &+phase_array(2)%volfracg, maxvolfrac
       if (I_AM_NODE_ZERO) then
          write (*,'(2(A25,2x,g17.8,/))')&
               'Re_m desired = ', Re, &
               'Re_m Output = ', test_mixmean_vmag*char_length/vis 
       endif
       
       test_mix_var = zero
       test_spec_var = zero
       pstart = 1
       do iphs = 1, nphases
          pend = pstart + phase_array(iphs)%npart - 1
          do m = pstart, pend
             do idim = 1, ndim
                test_spec_var(iphs) = test_spec_var(iphs) + (velbdy(m,idim)-test_mean_spec_vel(iphs,idim))**2.d0
             enddo
          enddo
          test_spec_var(iphs) = test_spec_var(iphs)/real(phase_array(iphs)%npart,prcn)
          pstart = pend + 1
       enddo
       do iphs = 1, nphases
          test_mix_var = test_mix_var + test_spec_var(iphs)*phase_array(iphs)%volfracg/maxvolfrac
       enddo
       if (I_AM_NODE_ZERO) then
          write (*,'(2(A25,2x,g17.8,/))')&
               'T desired = ', gran_temp, &
               'T Output = ',  test_mix_var/three
          
          write (*,'(2(A25,2x,g17.8,/))')&
               'Re_T desired = ', ReT, &
               'Re_T Output = ', DSQRT(test_mix_var/three)*char_length/vis 
       endif

       test_mix_var = zero
       test_spec_var = zero
       pstart = 1
       do iphs = 1, nphases
          pend = pstart + phase_array(iphs)%npart- 1
          do m = pstart, pend
             do idim = 1, ndim
                test_spec_var(iphs) = test_spec_var(iphs) + (velbdy(m,idim)-test_mixmean_vel(idim))**2.d0
             enddo
          enddo
          test_spec_var(iphs) = test_spec_var(iphs)/real(phase_array(iphs)%npart,prcn)
          pstart = pend + 1
       enddo
       do iphs = 1, nphases
          test_mix_var = test_mix_var + test_spec_var(iphs)*phase_array(iphs)%volfrac/mean_volfrac
       enddo
       
!!$       write (*,'(2(A25,2x,g17.8,/))')&
!!$            'Re_T desired = ', ReT, &
!!$            'Re_T Output HRENYA deF= ', DSQRT(test_mix_var/three)*dia_phys/vis 
    elseif (trim(psd_type).eq."csd") then
       test_mixmean_vel = zero
       do idim = 1, ndim
          test_mixmean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
       enddo
       
       test_mixmean_vmag = DSQRT(doT_PRODUCT(test_mixmean_vel(1:ndim),test_mixmean_vel(1:ndim)))
       if (I_AM_NODE_ZERO) then
          write (*,'(2(A25,2x,g17.8,/))')&
               'Re_m desired = ', Re, &
               'Re_m Output = ', test_mixmean_vmag*char_length/vis 
       endif
       test_mix_var = zero
       do m = 1, nbody
          do idim = 1, ndim
             test_mix_var = test_mix_var + (velbdy(m,idim)- test_mixmean_vel(idim))**2.d0
          enddo
       enddo
       test_mix_var = test_mix_var/real(nbody,prcn)
       if (I_AM_NODE_ZERO) then       
          write (*,'(2(A25,2x,g17.8,/))')&
               'Re_T desired = ', ReT, &
               'Re_T Output = ', DSQRT(test_mix_var/three)*char_length/vis 
       endif
    endif

    if (I_AM_NODE_ZERO)  call screen_separator(80,'I')
  end subroutine INPUT_CHECK
         
  
	subroutine alloc_mem
		use bcsetarrays
		use nlmainarrays
		use global_data 
		use field_tmp_arrays

		implicit none 
		integer :: rg, npt, i,j,k, iphs
!!$    revgroup = 2 
!!$    allocate(revp(revgroup))
!!$    allocate(pgrev  (revgroup))
!!$
!!$    revp(1:2) = nrpr
!!$    do rg = 1,revgroup
!!$       npt = revp(rg)
!!$       allocate(pgrev  (rg)%p(npt))
!!$    enddo

		allocate(gstencil  (order,order,order,3), vstencil(order,order,order,3), prsten(order,order,order))

		allocate(vsten(order,order,order,3), nlsten(order,order,order,3)&
			&,onlsten(order,order,order,3), dfsten(order,order,order,3),&
			& ppgrsten(order,order,order,3))

		allocate( gradphisten(order, order, order, 3, 3))
		allocate(sourcesten(order,order,order,3))
		!print*,'In alloc_mem_initflo ', size(gstencil,1), size(gstencil,2)
		!read (*,*)
		!nlarrays 

		!global_data

		if (igeometry==0) then
#if !PARALLEL
			allocate(u(local_no(1),local_no(2),local_no(3),ndim))
			allocate(p(local_no(1),local_no(2),local_no(3)))
			if (use_fes) allocate(pstar(local_no(1),local_no(2),local_no(3)))

			allocate( nl(local_no(1),local_no(2),local_no(3),ndim))
			allocate(onl(local_no(1),local_no(2),local_no(3),ndim))
			allocate( ff(local_no(1),local_no(2),local_no(3),ndim))

			allocate(wx(local_no(1)), wy(local_no(2)), wz(local_no(3)))
			!if (aliasflag.eq.1) then
			!	allocate(shifty(my),shiftz(mz),shiftyz(my,mz))
			!endif
			allocate(w2(local_no(1),local_no(2),local_no(3))) !, g(mx))
			allocate(uftmp(local_no(1),local_no(2),local_no(3)))

#else
			allocate(u(local_no(3),local_no(1),local_no(2),ndim))
			allocate(p(local_no(3),local_no(1),local_no(2)))
			if (use_fes) allocate(pstar(local_no(3),local_no(1),local_no(2)))

			allocate( nl(local_no(3),local_no(1),local_no(2),ndim))
			allocate(onl(local_no(3),local_no(1),local_no(2),ndim))
			allocate( ff(local_no(3),local_no(1),local_no(2),ndim))

			allocate(wx(local_no(3)), wy(local_no(1)), wz(local_no(2)))
			!if (aliasflag.eq.1) then
			!	allocate(shifty(my),shiftz(mz),shiftyz(my,mz))
			!endif
			allocate(w2(local_no(3),local_no(1),local_no(2))) !, g(mx))
			allocate(uftmp(local_no(3),local_no(1),local_no(2)))
#endif

			allocate(   fr(0:local_ni(1)+1, 0:local_ni(2)+1, 0:local_ni(3)+1, ndim))
			allocate(  ppr(0:local_ni(1)+1, 0:local_ni(2)+1, 0:local_ni(3)+1, ndim))
			allocate(diffn(0:local_ni(1)+1, 0:local_ni(2)+1, 0:local_ni(3)+1, ndim))

			allocate(pbc(0:local_ni(1)+1, 0:local_ni(2)+1, 0:local_ni(3)+1))
			allocate(ubc(-1:local_ni(1)+2, -1:local_ni(2)+2, -1:local_ni(3)+2, ndim))
			!extra buffers needed to store the velocity for interpolating the velocity at external points

			allocate(nlbc (0:local_ni(1)+1, 0:local_ni(2)+1, 0:local_ni(3)+1,ndim))
			allocate(onlbc(0:local_ni(1)+1, 0:local_ni(2)+1, 0:local_ni(3)+1,ndim))

			allocate(urtmp(local_ni(1),local_ni(2),local_ni(3)))

			u(:,:,:,:)   = czero
			p(:,:,:)     = czero
			nl(:,:,:,:)  = czero
			onl(:,:,:,:) = czero
			ff(:,:,:,:)  = czero

			fr(:,:,:,:)  =  zero
			ppr(:,:,:,:) =  zero
			pbc(:,:,:)   =  zero

			diffn(:,:,:,:) = zero 
			onlbc(:,:,:,:) = zero
			nlbc(:,:,:,:)  = zero
			ubc(:,:,:,:)   = zero

			urtmp(:,:,:) = zero
			uftmp(:,:,:) = zero

			ubcp=>ubc 
			nlbcp=>nlbc
			onlbcp =>onlbc
			pbcp => pbc

			allocate(fluid_atijk(0:local_ni(1)+1,0:local_ni(2)+1,0:local_ni(3)+1))
			allocate(GNACC_PART(local_ni(1),local_ni(2),local_ni(3),maxgn_partneigh+1))

			fluid_atijk = .true.

			!CONSTRUCTING VECTORS FOR DATA COMMUNICATION: LOGICAL FOR FLUID_ATIJK
			!XY-PLANE
			CREATE_XY_LSLICE( (local_ni(1)+2)*(local_ni(2)+2), xy_logical)
			COMMIT(xy_logical)
			!XZ-PLANE
			CREATE_XZ_LSLICE( local_ni(3)+2, local_ni(1)+2, (local_ni(1)+2)*(local_ni(2)+2), xz_logical)
			COMMIT(xz_logical)

			!CONSTRUCTING VECTORS FOR DATA COMMUNICATION: DOUBLE FOR ALL QUANTITIES EXCEPT VELOCITY
			!XY-PLANE
			CREATE_XY_RSLICE( (local_ni(1)+2)*(local_ni(2)+2), xy_double)
			COMMIT(xy_double)
			!XZ-PLANE
			CREATE_XZ_RSLICE( local_ni(3)+2, local_ni(1)+2, (local_ni(1)+2)*(local_ni(2)+2), xz_double)
			COMMIT(xz_double)

			!CONSTRUCTING VECTORS FOR DATA COMMUNICATION: FOR VELOCITY INCLUDING GOHST CELLS
			!XY-PLANE
			CREATE_XY_RSLICE( (local_ni(1)+4)*(local_ni(2)+4), xy_vel)
			COMMIT(xy_vel)
			!XZ-PLANE
			CREATE_XZ_RSLICE( local_ni(3)+4, local_ni(1)+4, (local_ni(1)+4)*(local_ni(2)+4), xz_vel)
			COMMIT(xz_vel)
		endif

		allocate(ferror_array(nerr_steps))
		do iphs = 1, nphases
			allocate(phase_array(iphs)%ferror_array(nerr_steps))
			if (imove.eq.1) allocate(phase_array(iphs)%grant_array(nerr_steps))
		enddo
	end subroutine alloc_mem

	subroutine set_interpolation_scheme(choice)
		use global_data, ONLY : scheme, interp_scheme, order,ob2l,ob2r&
				&,gstencil,intx_per, inty_per,intz_per, vstencil, vsten,&
				& prsten, ppgrsten, nlsten, onlsten, dfsten
		 implicit none 
		 integer, intENT(in) :: choice
		 integer :: order_orig
		 order_orig = order

		if (choice.EQ.1) then 
			interp_scheme = 'lpi'
			scheme = '4-order'
		elseif (choice.EQ.2) then 
			interp_scheme = 'lpi'
			scheme = '2-order'
		elseif (choice.EQ.3) then 
			interp_scheme = 'csi'
			scheme = '4-order'
		endif
		SELECT CASE(scheme)
		CASE("2-order")
			order = 2
		CASE("3-order")
			order = 3
		CASE("4-order")
			order = 4
		CASE("5-order")
			order = 5
		CASE("6-order")
			order = 6
		end SELECT

    if (allocated(gstencil).AND.order_orig.NE.order) then 
       deallocate(gstencil) 
       allocate(gstencil  (order,order,order,3))
    endif

    if (allocated(vstencil).AND.order_orig.NE.order) then 
       deallocate(vstencil) 
       allocate(vstencil  (order,order,order,3))
    endif

    if (allocated(prsten).AND.order_orig.NE.order) then 
       deallocate(prsten) 
       allocate(prsten  (order,order,order))
    endif

    if (allocated(vsten).AND.order_orig.NE.order) then 
       deallocate(vsten) 
       allocate(vsten  (order,order,order,3))
    endif

    if (allocated(nlsten).AND.order_orig.NE.order) then 
       deallocate(nlsten) 
       allocate(nlsten  (order,order,order,3))
    endif

    if (allocated(onlsten).AND.order_orig.NE.order) then 
       deallocate(onlsten) 
       allocate(onlsten  (order,order,order,3))
    endif

    if (allocated(dfsten).AND.order_orig.NE.order) then 
       deallocate(dfsten) 
       allocate(dfsten  (order,order,order,3))
    endif

    if (allocated(ppgrsten).AND.order_orig.NE.order) then 
       deallocate(ppgrsten) 
       allocate(ppgrsten  (order,order,order,3))
    endif
    if (allocated(gradphisten).AND.order_orig.NE.order) then 
       deallocate(gradphisten) 
       allocate(gradphisten  (order,order,order,3,3))
    endif
    intx_per = .false.
    !new if (xperiodic) intx_per = .true.
    intx_per = .true.
    inty_per = .true.
    intz_per = .true.
    ob2l = (order+1)/2
    ob2r = order/2 
  end subroutine set_interpolation_scheme

	subroutine generate_psd_config
		use functions
		implicit none
		real(prcn) :: davg,sigmad,drms,lbydrms,phi_calc,drms_calc,min_dia,max_dia, davg_calc,sigma_calc
		real(prcn), allocatable, dimension(:) :: FRANDN, DBDY
		real(prcn), allocatable, dimension(:) :: FTEMP,Wt
		real(prcn), allocatable, dimension(:) :: dia_inc,vfrac_inc, radii
		real(prcn), allocatable, dimension(:) :: moments,weights,abscissa

		real(prcn) :: fmin,fmax,hist(15), sum_volfrac
		real(prcn) :: vbox, x1, x2, sigd3, sigd2, conf
		integer :: SIZE_FRANDN, IDIM, NCOUNT, M, unit1,iphs, partstart, partend, sum_part
		integer :: nsim, imis_gof, j

		character*100 filename
		real(prcn) :: confint, int_dist_avg, int_dist_var, int_dist_sd, int_dist2


		davg = dia_phys
		sigmad = sigmabydavg*davg

		if (trim(psd_type).eq.'csd') nphases = UNdeFINED_I
		if (nphases.eq.1) sigmad = zero

		drms = davg*(DSQRT(1 + (sigmad*davg)**2))

		!lbydrms = doml(2)/drms
		!print*,'DRMS = ', DRMS

		if (I_AM_NODE_ZERO) then
			if (trim(psd_type).eq.'csd') then
				char_length = davg
				doml(2) = lybyd*davg
				doml(1) = doml(2)
				doml(3) = doml(2)

				print*,'LYBYD = ', LYBYD, 'DAVG = ', davg, 'doml(2) = ', doml(2)
				vbox = doml(1)*doml(2)*doml(3)

				nbody = Nint(phiavg*6.d0*(lybyd**3.d0)/pi)
				if (I_AM_NODE_ZERO) print*,'NBODY IN CSD CASE  = ', NBODY
				SIZE_FRANDN = 4*NBODY
				allocate(FRANDN(SIZE_FRANDN),DBDY(nbody))
599			CONTINUE
				call norm_dist(FRANDN(1:SIZE_FRANDN))
				ncount = 0 
				do m = 1, SIZE_FRANDN
					if (ABS(FRANDN(m)).lt.3.d0) then
						ncount = ncount + 1
						DBDY(NCOUNT) = FRANDN(M)*sigmad + davg
					endif
					if (NCOUNT.eq.NBODY)EXIT
				enddo
				if (NCOUNT.lt.NBODY) then
					write (*,*)'NCOUNT = ', 'ncount, NBODY = ', NBODY
					write (*,*)'REdoING NORMAL DISTRIBUTION'
					GOTO 599
				endif
				min_dia = MINVAL(DBDY)
				MAX_DIA = MAXVAL(DBDY)
				DAVG_CALC = ZERO
				DRMS_CALC = ZERO
				PHI_CALC = ZERO
				do m = 1, nbody
					davg_calc = davg_calc + dbdy(m)/real(nbody,prcn)
					drms_calc = drms_calc+ (dbdy(m))**2
					phi_calc = phi_calc + pi*dbdy(m)*(dbdy(m))**2/6.d0
				enddo
				drms_calc = DSQRT(drms_calc/real(nbody,prcn))
				phi_calc = phi_calc/(doml(1)*doml(2)*doml(3))
				sigma_calc = DSQRT((drms_calc)**2- (davg_calc)**2)
				if (ABS(sigma_calc-sigmabydavg*davg).gt.0.01*sigmabydavg*davg) then
					write (*,'(A80)')'REQUIRED SIGMA/<D> NOT REACHED. doING NORM DIST AGAIN'
					goto 599
				endif
				write (*,'(3(2x,A25,g17.8,/))')'MIN DIA = ', min_dia, 'MAX DIA =', max_dia,'RATIO = ', max_dia/min_dia
				write (*,'(2(2x,A25,g17.8,/))')'DIA AVG CALC = ', davg_calc, 'DIA AVG INPUT =', davg

				write (*,'(2(2x,A25,g17.8,/))')'VOL FRAC CALC = ', phi_calc, 'VOL FRAC INPUT =', phiavg
				write (*,'(2(2x,A25,g17.8,/)))')'DRMS CALC = ', drms_calc, 'DRMS INPUT =', drms
				write (*,'(2(2x,A25,g17.8,/)))')'SIGMA/<D> CALC = ', DSQRT((drms_calc/davg_calc)**2-one), 'SIGMA/<D> INPUT =', sigmabydavg
				!!$ nit1=getnewunit(minunitno,maxunitno)
				!!$ pen (unit=unit1,file='psd.dat',status='unknown')
				!!$ allocate(ftemp(nbody),wt(nbody))
				!!$ do m=1,nbody
				!!$ wt(m) = 1/real(nbody,prcn)
				!!$ enddo
				!!$ ftemp(1:nbody) = dbdy(1:nbody)
				!!$ call histogram(ftemp,wt,nbody,15,fmin,fmax,hist)
				!!$ call plothist(hist(1:15),fmin,fmax,15,unit1,1.d0,1.d0)
				!!$ close(unit1,status='keep')
				nphases = nbody
				call alloc_phase_related
				do iphs = 1, nphases
					phase_array(iphs)%npart = 1
					phase_array(iphs)%dia = dbdy(iphs)
					phase_array(iphs)%volfrac = (phase_array(iphs)%npart)*pi*(phase_array(iphs)%dia**3.d0)/(6.d0*vbox)
				enddo
			elseif (trim(psd_type).eq.'discrete') then
				write (*,'(A,2x,I6)')'GENERATING A REPRESENTATION OF GAUSSIAN CSD WITH PHASES = ', nphases
				allocate(moments(2*nphases),weights(nphases),abscissa(nphases))

				char_length = davg
				doml(2) = lybyd*davg
				doml(1) = doml(2)
				doml(3) = doml(2)

				print*,'LYBYD = ', LYBYD, 'DAVG = ', davg, 'doml(2) = ', doml(2)
				vbox = doml(1)*doml(2)*doml(3)

				nbody = int(phiavg*6.d0*(vbox/davg**3.d0)/pi)
				char_length = davg
				if (I_AM_NODE_ZERO) print*,'NBODY IN DISCRETE CASE  = ', NBODY
				moments(:) = zero
				moments(1) = 1
				if (nphases.gt.1) then
					moments(3) = 1
					do iphs = 3,nphases
						moments(2*iphs-1) = moments(2*iphs-3)*(2*iphs-3)
					enddo
				endif
				print*,'MOMENTS =', moments(:)
				call hermite_polynomial(nphases,moments,weights,abscissa) 
				call alloc_phase_related
				allocate(DBDY(nbody))
				sum_part = 0
				do iphs = 2, nphases
					phase_array(iphs)%npart = int(weights(iphs)*nbody)
					sum_part = sum_part+phase_array(iphs)%npart
				enddo
!!$			npart(nphases) = nbody - SUM(npart(1:nphases-1))
				if (nphases.gt.1) then
					phase_array(1)%npart = nbody - sum_part
				else
					phase_array(1)%npart = nbody
				endif
				sum_volfrac = zero
				do iphs = 1, nphases
					phase_array(iphs)%dia = davg + abscissa(iphs)*sigmad
					phase_array(iphs)%volfrac = phase_array(iphs)%npart*pi*(phase_array(iphs)%dia**3.d0)/(6.d0*vbox)
					sum_volfrac = sum_volfrac + phase_array(iphs)%volfrac
				enddo
          
				write (*,'(2(2x,A25,g17.8,/))')'VOL FRAC CALC = ', sum_volfrac, 'VOL FRAC INPUT =', phiavg
				write (*,'(2(2x,A25,I6,/))') 'NBODY CALC = ', sum_part+phase_array(1)%npart, 'NBODY INPUT =', nbody
			elseif (trim(psd_type).eq.'bidisp') then
				nphases = 2
				ibidisperse = .true.
				call alloc_phase_related
				dia_ratio = yalpha(2)/yalpha(1)
				volfrac_rat = (dia_ratio-yalpha(2))/(yalpha(2)-1)

				write (*,'(A)')'GENERATING A BIDISPERSE SUSPENSION WITH : '
				write (*,'(A, 2x, g17.8)')'VOLUME FRACTION RATIO = ', volfrac_rat
				write (*,'(A, 2x, g17.8)')'DIAMETER RATIO = ', dia_ratio

				phase_array(1)%volfrac = phiavg/(one+volfrac_rat)
				phase_array(2)%volfrac = volfrac_rat*(phase_array(1)%volfrac)

				x1 = phase_array(1)%volfrac/phiavg
				x2 = phase_array(2)%volfrac/phiavg

				phase_array(1)%dia = dia_phys*(x1 + x2/dia_ratio)
				phase_array(2)%dia = dia_ratio*(phase_array(1)%dia)

				!doml(2) = lybyd*phase_array(2)%dia
				doml(2) = lybyd*dia_phys
				doml(1) = doml(2)
				doml(3) = doml(2)
				print*,'LYBYD = ', LYBYD, 'DIA2 = ',phase_array(2)%dia , 'doml(2) = ', doml(2)
				vbox = doml(1)*doml(2)*doml(3)

				phase_array(1)%npart = int(6.d0*vbox*phase_array(1)%volfrac/(pi*phase_array(1)%dia**3.d0))
				phase_array(2)%npart = int(6.d0*vbox*phase_array(2)%volfrac/(pi*phase_array(2)%dia**3.d0))

				nbody = phase_array(1)%npart + phase_array(2)%npart

				!^^^^^^SINGLE PARTiCL MOVING ^^^^^^^^
				!phase_array(1)%dia = dia_phys/3
				!phase_array(2)%dia = dia_phys
				!phase_array(1)%volfrac = pi/6 * (phase_array(1)%dia/doml(2))**3
				!phase_array(2)%volfrac = phiavg
				!phase_array(1)%npart = 1
				!phase_array(2)%npart = int(6.d0*vbox*phase_array(2)%volfrac/(pi*phase_array(2)%dia**3.d0))
				!nbody = phase_array(1)%npart + phase_array(2)%npart
				!phiavg = phase_array(1)%volfrac + phase_array(2)%volfrac
				!------------------------------------
          
				write (*, '(A25,2(2x,g17.8))')'Diameters (ORIGINAL)=', phase_array(1)%dia, phase_array(2)%dia
				write (*, '(A25,2(2x,g17.8))')'Volume fractions (ORIGINAL)=', phase_array(1)%volfrac, phase_array(2)%volfrac
				write (*, '(A25,2(2x,I6))') 'No. of Particles =', phase_array(1)%npart, phase_array(2)%npart
				write (*,'(2(2x,A25,I6,/))') 'NBODY CALC = ', phase_array(1)%npart+phase_array(2)%npart

				allocate(DBDY(nbody))
			endif
			allocate(dia_inc(nphases),vfrac_inc(nphases), xc_gener(nbody,ndim), rad_gener(nbody))

			partstart = 1
			sum_volfrac = zero
			min_dia = LARGE_NUMBER
			do iphs = 1, nphases
				if (phase_array(iphs)%dia.lt.min_dia) then
					min_dia = phase_array(iphs)%dia
				endif
			enddo

			if (my.eq.undefined_I) then
				if (trim(psd_type).eq.'bidisp') then
					my = (doml(2)/dia_phys)*dbydx
				else
					my = (doml(2)/min_dia)*dbydx
				endif
				print*,'SETTING MY BASED ON THE dbydx FOR THE SMALLEST SPHERE, WHICH IS:', MY
				if (MOD(my,2).ne.0) then
					my = my+1
					print*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
				endif
				mx = my !new +1
				mz = my
				!mxf = mx
				!my2 = my/2+1
				!mx1 = mx !new - 1 
				!mx2 = mx/2+1
			endif
			dy = doml(2)/real(my,prcn)
			dx = dy
			dz = dy 

			do iphs = 1,nphases
				percent_buf(iphs) = (min_part_sep*dx/phase_array(iphs)%dia)*100.d0
				percent_buf(iphs) = percent_buf(iphs) + one

				dia_inc(iphs) = phase_array(iphs)%dia*(one + percent_buf(iphs)/100.d0)
				partend = partstart + phase_array(iphs)%npart-1
				do m = partstart,partend
					dbdy(m) = dia_inc(iphs)
				enddo
				partstart = partend + 1
				sum_volfrac = sum_volfrac + phase_array(iphs)%volfrac
				vfrac_inc(iphs) = phase_array(iphs)%volfrac*(one+percent_buf(iphs)/100.d0)**3.d0
				!print*,'percent_buf IN = ', percent_buf(iphs)
			enddo
       
			write (*,'(3(2x,A25,g17.8,/))')&
			'VOL ORIG = ', sum_volfrac, &
			'VOL FRAC INC =', SUM(vfrac_inc(:))

			if (GOF_AVG) then 
				write (*,*)'GOF AVG IS true, THEREFORE GENERATING', mis_gof, 'initializations, averaging and then stopping' 
				write (ounit,*)'GOF AVG IS true, THEREFORE GENERATING', mis_gof, 'initializations, averaging and then stopping' 
				nsim = mis_hs
			else 
				nsim = 1
			endif

			if (.not.allocated(gofr_avg)) then
				allocate(gofr_avg(nbins), gofr_mis(nsim, nbins))
				allocate(rad_bin(nbins))
				allocate(int_dist(nsim))

				gofr_avg = zero
				gofr_mis = zero
				rad_bin  = zero
				int_dist = zero
			endif

			do imis_gof = 1, nsim
				write (*,*) 'GENERATING RANdoM CONFIG FOR IMIS =', imis_gof		 
				call gener_config_random(dbdy(1:nbody),dia_inc(1:nphases),vfrac_inc(1:nphases))
				! gofr_mis(imis_gof,:) = gofr(:)
				! call interstitial_dist(xc_gener, rad_gener, int_dist(imis_gof), nbody)
			enddo

#if 0
			gofr_avg = zero 
			do j=1,nbins
				gofr_avg(j) = gofr_avg(j) + one/dble(nsim)*sum(gofr_mis(1:nsim,j))
			enddo

			do j=1,nbins
				if (nsim.ge.2) then 
					conf=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* sum((gofr_mis(1:nsim,j)- gofr_avg(j))**2))
				endif

				write (gofavgunit,'(10(E20.10,1x))') rad_bin(j), rad_bin(j)*lybyd, gofr_avg(j), conf
			enddo
			close(gofavgunit, status = "keep") 

			int_dist_avg = sum(int_dist)/nsim
			if (nsim>2) then
				call get_confin(nsim, confint)
				int_dist_var = sum( (int_dist(1:nsim)- int_dist_avg) **2) / nsim
				int_dist_sd  = sqrt(int_dist_var)
				conf = int_dist_sd / sqrt(float(nsim)) * confint
			else
				conf = zero
			endif

			!call int_dist_gofr(nbins, rad_bin(1:nbins), gofr_avg(1:nbins), int_dist2)
			filename = "NMIS_interparticle.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")
			write (1,"(4D15.7)") maxvolfrac, int_dist_avg, conf, int_dist2
			close (1)
#endif

			if (GOF_AVG) then 
				if (I_AM_NODE_ZERO) write (*,*) 'GOF AVG IS true, SO stopPING THE SIMULATION'
				PARALLEL_FINISH()
				stop 
			endif

			if (trim(psd_type).eq.'csd') then
				nphases = nbody
				call dealloc_phase_related
				call alloc_phase_related
				do iphs = 1, nphases
					phase_array(iphs)%npart = 1
				enddo
			endif
		endif
    
		BROADCAST_INT(nphases,1,node_zero,comm_cart_2d)

		if (.not.I_AM_NODE_ZERO) call alloc_phase_related
		sum_part = 0
		do iphs = 1, nphases
			BROADCAST_INT(phase_array(iphs)%npart,1,node_zero,comm_cart_2d)
			sum_part = sum_part + phase_array(iphs)%npart
		enddo
		nbody = sum_part

		call alloc_bnd_related
		if (I_AM_NODE_ZERO) then
			do m = 1, NBODY
				xc(m,1:ndim) = xc_gener(m,1:ndim)
				RADBDY(m) = RAD_GENER(m)
			enddo
			deallocate(XC_GENER, RAD_GENER)
		endif

		BROADCAST_INT(mx,1,node_zero,comm_cart_2d)
		BROADCAST_INT(my,1,node_zero,comm_cart_2d)
		BROADCAST_INT(mz,1,node_zero,comm_cart_2d)
		BROADCAST_DOUBLE(doml(1),ndim,node_zero,comm_cart_2d)

		if (.not.I_AM_NODE_ZERO) then
			dy = doml(2)/real(my,prcn)
			dx = dy
			dz = dy 
			!mx = my !new +1
			!mz = my
			!mxf = mx
			!my2 =  my/2+1
			!mx1 = mx !new - 1 
			!mx2 = mx/2+1
		endif

		if (I_AM_NODE_ZERO) then
			write (*,*) 'BOX SIZE   = ', mx, my, mz
			!write (*,*) 'mxf, mx1, my2   = ', mxf, mx1, my2
			if (trim(psd_type).eq.'bidisp') then
				sigd3 = zero
				sigd2 = zero
				do m = 1, nbody
					sigd3 = sigd3 + (two*radbdy(m)*dx)**3.d0
					sigd2 = sigd2 + (two*radbdy(m)*dx)**2.d0
				enddo
				char_length = sigd3/sigd2
				write (*,'(A30, 2x,g17.8)')'CHAR LENGTH FOR BIDISP CASE = ', char_length

				if (I_AM_NODE_ZERO) then
					open (unit = 1001, file = trim(run_name)//'_yis.dat', status='unknown')
					write (1001,*)phase_array(1)%dia/char_length, phase_array(2)%dia/char_length
					close(1001, status='keep')
				endif
			endif
		endif 
#if PARALLEL
		BROADCAST_DOUBLE(char_length, 1, node_zero, comm_cart_2d)
#endif
	end subroutine generate_psd_config


	subroutine gener_config_random(DBDY,DIA_IN, PHI_IN)
		use dem_mod
		use collision_mod
		use dependent_functions
		implicit none 

		real(prcn), intENT(IN),dimension(:) :: DIA_IN, PHI_IN,DBDY
		real(prcn) :: vbox, xtemp, ly_lat, dmax, davg,&
				& DAVG_CALC, SIGMAD_CALC, Tstop1, Tstop2, YLEN,&
				& doml_TMP(3), A, B, ymin, ymax , DIST&
				&, davgmax, sigmabydavgmax,  DAVG_BAR, RHOP_BAR, TMP&
				&,dia_sm , DIA_AIM, RHO_AIM, DRATIO, mfp, tmfp
		real(prcn) :: total_phi, min_dia, max_dia, drms_calc, ndens
		real(prcn), dimension(:), allocatable :: frandn, radbdy_tmp
		integer ::  nratio, MCOUNT, MINdeX, NSPEC_AIM,&
				& i,j,k, NP, iphs, partstart, partend
		integer SIZE_FRANDN, IDIM, NSPEC_INIT(nphases), trialcount, npart(nphases)
		character*80 ::  collision_type_orig    
		integer, dimension(1):: mloc
		logical, allocatable :: contact(:,:)
		real(prcn) :: max_overlap


		call screen_separator(80,'R')
		vbox = doml(1)*doml(2)*doml(3)       
		collision_type_orig = collision_type
		do iphs = 1, nphases
			npart(iphs) = phase_array(iphs)%npart
		enddo

		nspec_init(1:nphases) = npart(1:nphases)

		SIZE_FRANDN = 4*NBODY

		allocate(radbdy_TMP(nbody), frandn(SIZE_FRANDN))

		DMAX = MAXVAL(DIA_IN(1:nphases))

		trialcount = 1
!!$    
!!$    if (NSPEC1.LT.NSPEC2) then
!!$       dbdy(:) = dia2_in
!!$       NSPEC_AIM = NSPEC1
!!$       DIA_AIM = dia1_in
!!$    else
!!$       dbdy(:) = dia1_in
!!$       NSPEC_AIM = NSPEC2
!!$       DIA_AIM = dia2_in
!!$    endif
!!$    
!!$    if (.not.ibidisperse) goto 2000
!!$    
!!$1000 CONTINUE
!!$    !       if (DRATIO.GT.ONE) then 
!!$    MCOUNT = 0
!!$    call uni_dist(FRANDN(1:SIZE(FRANDN)))
!!$    FRANDN(:) = FRANDN(:)*NBODY
!!$
!!$    do I = 1, SIZE(FRANDN)
!!$       
!!$       MINdeX = Nint(FRANDN(I))!NBODY
!!$       if (MINdeX.EQ.0) CYCLE
!!$!       print*,'MINdeX = ', MINdeX, DIA_AIM, NSPEC_AIM
!!$       if (DRATIO.GT.ONE) then 
!!$          if (dbdy(MINdeX).NE.DIA_AIM) then 
!!$             DBDY(MINdeX) = DIA_AIM
!!$             MCOUNT = MCOUNT + 1
!!$             if (MCOUNT.EQ.NSPEC_AIM) GOTO 2000
!!$          endif
!!$       else
!!$          MCOUNT = MCOUNT + 1
!!$          if (MCOUNT.EQ.NSPEC_AIM) GOTO 2000
!!$       endif
!!$    enddo
!!$    
!!$    if (MCOUNT.LT.NSPEC_AIM) then 
!!$       trialcount = trialcount+1
!!$       print*,'doING UNIDIST AGAIN'
!!$       print*,'NSPEC_AIM = ', NSPEC_AIM, ' AND MCOUNT = ', MCOUNT
!!$       if (trialcount.gt.50) then 
!!$          write (*,'(A,/,A,/,A)')'EXCEEdeD 50 TRIALS OF TRYING TO GENERATE MINdeX',' GIVING UP AND stopPING HERE', ' TRY A DifFERENT INITIAL SEED.'
!!$          stop
!!$       else
!!$          GOTO 1000
!!$       endif
!!$    endif
!!$          
!!$2000 CONTINUE
!!$    !       endif
!!$    DAVG = ZERO
!!$    do I = 1, NBODY
!!$       DAVG  = DAVG + DBDY(I)/NBODY
!!$    enddo
!!$

!!$
!!$    do I = 1, NBODY
!!$       TOTAL_PHI = TOTAL_PHI + PI*(DBDY(I)**THREE)/6.d0
!!$       DRMS_CALC  = DRMS_CALC + DBDY(I)**2.d0
!!$       DAVG_CALC  = DAVG_CALC + DBDY(I) 
!!$       
!!$    enddo
!!$
!!$    
!!$    TOTAL_PHI = TOTAL_PHI/(doml(1)*doml(2)*doml(3))
!!$
!!$    print*,'TOTAL VOL FRACTION', TOTAL_PHI
!!$    write (*,*) 'MAX DIA RATIO = ', MAXVAL(DBDY)/MINVAL((DBDY))
    
		doml_TMP(:) = doml(:)
		MAX_DIA = MAXVAL(DBDY)
		min_dia = MINVAL(DBDY)
		print*,'DMAX = ', MAX_DIA, min_dia
		call gener_lattice_mod(nbody,doml(1:3),xc_gener(1:nbody,1:ndim),dmax,dbdy(1:nbody))
		print*,'OUT OF GENER_LATTICE_MOD'
		if (MAXVAL(doml(1)-XC_gener(1:nbody,1))-MAX_DIA*HALF.LT.ZERO) print*,'PARTICLE COMING OUT OF +X'
		if (MINVAL(XC_gener(1:nbody,1))-MAX_DIA*HALF.LT.ZERO) then 
			print*,'PARTICLE COMING OUT OF -X', MINVAL(XC_gener(1:nbody,1)), MAX_DIA, MINVAL(XC_gener(:,1))- MAX_DIA
		endif
    
		if ((doml(2)-MAXVAL(XC_GENER(1:nbody,2)))-MAX_DIA*HALF.LT.ZERO) print*,'PARTICLE COMING OUT OF +Y'
		if (MINVAL(XC_GENER(1:nbody,2))-MAX_DIA*HALF.LT.ZERO) print*,'PARTICLE COMING OUT OF -Y'
		if (MAXVAL(doml(3)-XC_GENER(1:nbody,3))-MAX_DIA*HALF.LT.ZERO) print*,'PARTICLE COMING OUT OF +Z'
		if (MINVAL(XC_GENER(1:nbody,3))-MAX_DIA*HALF.LT.ZERO) print*,'PARTICLE COMING OUT OF -Z'
    
		XLENGTH = doml(1)
		ZLENGTH = doml(3)
		YLENGTH  =  MAXVAL(XC_GENER(1:nbody,2))+DMAX*HALF

   
		RAD_GENER(1:NbodY) = HALF*DBDY(1:NBODY)
		!goto 20000
		GENER_CONFIG_CASE = .true.    
		DT = LARGE_NUMBER
    
		if (YLENGTH.LT.doml(2)) then 
			write (*,'(A7,2x,g17.8,A7,2x,g17.8)')'YLEN  = ', YLENGTH, 'LT doml(2) = ', doml(2)

			SHRINK = .false.

			ymin = MINVAL(XC_GENER(1:nbody,2)) - MINVAL(dbdy(1:nbody))*half
			ymax = MAXVAL(XC_GENER(1:nbody,2)) + MAXVAL(dbdy(1:nbody))*half

			A = doml(2)/(ymax-ymin)
			B = -A*ymin
			print*,'A, B = ', A,b
			XC_GENER(1:nbody,2) = XC_GENER(1:nbody,2)*A + B
			!call des_time_march(.true.)    
		else 
			write (*,'(A7,2x,g17.8,A7,2x,g17.8)')'YLEN  = ', YLENGTH, 'GT doml(2) = ', doml(2)
			collision_type="softsphere"
			Tstop = 20.d0*SQRT((two*YLENGTH)/980.d0)
			SHRINK = .true.
			TEST_YMAXVAL = .true.
			YMAXVAL = doml(2)
			deS_EN_INPUT = 0.3
			deS_EN_WALL_INPUT = 1.0

			call  des_time_march(.true.)
			call  des_time_march(.false.)

			SHRINK = .false.
			TEST_YMAXVAL = .false.

			XC_GENER(1:NBODY,1:3) = deS_POS_NEW(1:NBODY, 1:3)
		endif
    
		collision_type="softsphere"
    
		XLENGTH = doml(1)
		ZLENGTH = doml(3)
		YLENGTH = doml(2)
		deS_EN_INPUT(:) = 1.0
		deS_EN_WALL_INPUT(:) = 1.0
		!print*,'SI =', SIZE(deS_EN_INPUT), deS_EN_INPUT
		!!$    ndens = (nspec1+nspec2)/(vbox)
		ndens = SUM(npart(1:nphases))/(vbox)
		!print*,'NUMBER deNSITY = ', ndens
		mfp = (one/ndens)**(one/three)
		!mfp = dchar/(phi1+phi2)
		tmfp = mfp/dsqrt(pvel_var)
		tstop = 10.d0*tmfp

		write (*,*)'mfp, tmfp, tstop = ', mfp, tmfp, tstop
!20000 continue
    
		call  des_time_march(.true.)
		call  des_time_march(.false.)

		do IDIM = 1, 3
			deS_POS_NEW(1:NBODY,IDIM) = deS_POS_NEW(1:NBODY, IDIM)/doml(IDIM)
		enddo
    
		radbdy_tmp(1:NBODY)=rad_gener(1:NBODY)/doml(2)
		do I = 1, NBODY
			mloc = MINLOC(radbdy_tmp)
			rad_gener(I) = radbdy_tmp(MLOC(1))
			xc_gener(I,:) = des_pos_new(MLOC(1),:)
			radbdy_tmp(mloc(1)) = LARGE_NUMBER
		enddo
    
		call scale_to_grid_units(nbody,npart(1:nphases),nphases,my,my,xperiodic(1)&
			&,percent_buf(1:nphases),xc_gener(1:nbody,1:3),rad_gener(1:nbody),&
			& min_part_sep, toscale=.true.) 

		nbody = SUM(npart(1:nphases))
    
		do iphs = 1, nphases
			phase_array(iphs)%npart = npart(iphs)
			if (trim(psd_type).ne.'csd') then
				write (*,'(A80,g17.8,2x,i4)')'FINAL VOL FRAC AND NSPECS AFTER LATTICE, SHRINKAGE, AND RESCALING = ', npart(iphs)*fourthirdpi*(phase_array(iphs)%dia/(two*dx*real(my,prcn)))**3, npart(iphs)
				!!$ write (*,'(A80,g17.8,2x,i4)')'RADBDY = ', radbdy(1)
				write (*,'(A80,2(2x,i2))') 'NUMBER OF PARTICLES LOST FOR PAHSE = ', iphs, nspec_init(iphs)-npart(iphs)
				!read (*,*)
				!!$    if (nspec_init(iphs)-npart(iphs)1.gt.0.OR.nspec2_init-nspec2.gt.0) then 
				!!$	open (1001, file=trim(run_name)//"_deLETION_true.dat", form="formatted")
				!!$    write (*,'(A80,2(2x,i2))') 'NUMBER OF PARTICLES LOST IN 1st AND 2nd SPECIES = ', nspec1_init-nspec1, nspec2_init-nspec2
				!!$    close(1001, status="keep")
				!!$    endif
			endif
		enddo
    
!		nbins = 200
    
		rescaling = .true.
    
		open (1001, file=trim(run_name)//"_xc_post_grid_scal.dat", form="formatted")
		write (1001,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" ',' "MARK" '
		partstart = 1
		do iphs = 1, nphases
			partend = partstart + npart(iphs) -1 
			do NP = PARTSTART,PARTend
				write (1001,'(10(2x,g15.8))')( XC_GENER(NP, i), i = 1, ndim), rad_gener(NP),iphs
			enddo
			partstart = partend + 1
		enddo

		close(1001,status="keep")
    
#if 0
		if (xperiodic) then
			if (.not.allocated(contact)) allocate(contact(nbody,nbody))
			call calculate_gofr_homog(nbody,xc_gener(1:nbody,1:3), contact, global_n(:), nbins, rescaling, gofr(1:nbins), rad_bin(1:nbins), max_overlap)
		else
			call calc_gofr(nbody, xc_gener(1:nbody,1:3), rad_gener(1:nbody), my, mxf,xperiodic,nbins, 3, rescaling, gofr(1:nbins),rho_est(1:nbins), rad_bin(1:nbins)) 
		endif
#endif
    
		!!$    do j=1,nbins
		!!$       write (gofunit,'(10(E20.10,1x))')rad_bin(j),rad_bin(j)/(one/(two*lybyd))&
		!!$            &,rho_est(j),gofr(j)
		!!$    enddo
		!!$    close(gofunit, status = "keep")
    
		GENER_CONFIG_CASE = .false.
		collision_type = collision_type_orig
		call screen_separator(80,'R')
	end subroutine gener_config_random



  subroutine jacinv(nx,ny,nz,jj,cp,sp,st,ct,r1,r2)

    implicit none

    real(prcn) ::  nx,ny,nz
    integer ::  i,j
    real(prcn) ::  jj(3,3)
    real(prcn) :: st,ct,sp,cp

    real(prcn) :: r1,r2


    !	r1=dsqrt(ny*ny+nz*nz)
    !	r2=dsqrt(nx*nx+ny*ny+nz*nz)

    !	ct=ny/r1
    !	st=nz/r1

    if (r1.NE.(0.d0)) then

       !       ct=ny/r1
       !       st=nz/r1

       !       cp=nx/r2
       !       sp=r1/r2
       jj(1,1)=cp
       jj(1,2)=sp*ct
       jj(1,3)=sp*st

       jj(2,1)=-sp
       jj(2,2)=cp*ct
       jj(2,3)=cp*st

       jj(3,1)=0.d0
       jj(3,2)=-st/sp
       jj(3,3)=ct/sp

    else

       !special case for phi=zero
       !otherwise matrix is singular

       do i=1,3
          do j=1,3
             jj(i,j)=0.d0
          enddo
       enddo

       jj(1,1)=nx/ABS(nx)
       !write (*,*) nx,jj(1,1)
       jj(2,2)=1.d0
       jj(3,3)=1.d0

    endif

    return

  end subroutine jacinv




  
  subroutine alloc_bnd_related
    implicit none 
    
    integer :: NBODY_TMP

    allocate(force(nbody,ndim),pres(nbody,ndim),visc(nbody,ndim),torq(nbody,ndim), force_loc(nbody,ndim), force_chem(nbody,ndim), contact_force(nbody,ndim))
    allocate(presloc(nbody,ndim), viscloc(nbody,ndim),torqloc(nbody,ndim))

    allocate(ap(nbody,ndim,2),up(nbody,ndim,2),angv(nbody,ndim,2),anga(nbody,ndim,2))

    allocate(mp(nbody),mpart(nbody),mompart(nbody), xc(nbody, ndim), tmpfor(nbody,ndim), velbdy(nbody,ndim))

	allocate(color(nbody))

	angv = zero
	anga = zero

    allocate(radbdy(nbody),radibdy(nbody),radobdy(nbody),rado2bdy(nbody))
    allocate(part_array(nbody))
#if 0
    if (iscalon.eq.1) then
       allocate(phisurfall(nbody, nspmx))
       allocate(phirmean(nspmx), fphirmean(nspmx),sum_flux_nm1(nbody)&
            &, sum_flux_nm2(nbody), sum_flux(nbody), flux_body(nbody&
            &,nspmx), flux_body2(nbody,nspmx))
       phisurfall = zero
       phisurfall_nm1 = zero
       phisurfall_nm2 = zero
       sum_flux = zero
       sum_flux_nm1 = zero
       sum_flux_nm2 = zero
    endif
#endif
    
  end subroutine alloc_bnd_related

	subroutine dealloc_bnd_related
		implicit none 
		integer :: m

		deallocate(force,pres,visc,torq,force_loc,force_chem,contact_force)
		deallocate(presloc, viscloc,torqloc)
		deallocate(ap,up,angv,anga)
		deallocate(mp,mpart,mompart, xc, tmpfor, velbdy)

		deallocate(radbdy,radibdy,radobdy,rado2bdy)
		do m = 1, nbody
			if (associated(part_array(m)%if_rev))  deallocate(part_array(m)%if_rev)
			if (associated(part_array(m)%if_drag)) deallocate(part_array(m)%if_drag)
		enddo
		deallocate(part_array)

#if 0
		if (iscalon.eq.1) then
			deallocate(phisurfall)
			deallocate(phirmean, fphirmean,sum_flux_nm1, sum_flux_nm2, sum_flux, flux_body, flux_body2)
		endif
#endif
	end subroutine dealloc_bnd_related


  
  subroutine alloc_phase_related
    implicit none
    allocate(phase_array(nphases), percent_buf(nphases),norm_drag_spec(nphases),norm_drag_chem_spec(nphases))
  end subroutine alloc_phase_related

  subroutine dealloc_phase_related
    implicit none
    integer :: iphs

    do iphs = 1, nphases
       if (associated(phase_array(iphs)%bndpts)) deallocate(phase_array(iphs)%bndpts)
       if (associated(phase_array(iphs)%ferror_array)) deallocate(phase_array(iphs)%ferror_array)
    enddo

    deallocate(phase_array, percent_buf, norm_drag_spec, norm_drag_chem_spec)
  end subroutine dealloc_phase_related


#if 0
  subroutine doMAIN_1D_deCOMP( n, numprocs, myid, s, e )
    integer,intent(in):: n, numprocs, myid
    integer, intent(out) :: s, e
    integer :: nlocal, deficit
    
    nlocal  = n / numprocs
    s	      = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s	      = s + min(myid,deficit)
    
    if (myid .lt. deficit) then
       nlocal = nlocal + 1
    endif
    
    e = s + nlocal - 1
    
    if (e .gt. n .or. myid .eq. numprocs-1) e = n
    
  end subroutine doMAIN_1D_deCOMP
#endif

end MODULE initialize_flo
