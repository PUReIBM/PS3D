MODULE boundary_condition
#include "ibm.h"
	USE precision 
	use constants 
	USE global_data
	!new USE fftw_interface
	USE fftw3_interface      
	use parallel
	USE dependent_functions
	!USE  initialize_flo
	USE collision_mod
	!uncomment the following line omega needs to be calculated 
	Use nlmainarrays, Only : ur=>ubcp,pr=> pbcp!, nlr=>nlbcp, onlr=>onlbcp 
	USE bcsetarrays, ONLY :  fr, ppr, diffn 
	!USE nlarrays, ONLY : ff1=>uf1,ff2=>uf2,ff3=>uf3, fr1=>ur1, fr2=>ur2, fr3=>ur3, dudy=>ur11, dudz=>ur22, dvdz=>ur33, dwdy=>ur12
	use report
	use field_tmp_arrays
	use init_turb, only : forced_turbulence, linear_forcing, eddy_time_i, u_prime, activate_forcing_time
	implicit none 

	PRIVATE 
	real(prcn) ::  vort(ndim), da(2)
	real(prcn) ::  usph(ndim),ucar(ndim) 
	real(prcn) ::   furf

	real(prcn) ::  pl,nll(ndim),onll(ndim),ppll(ndim),dfll(ndim)
	real(prcn) ::  snorm(ndim),rad
	real(prcn) ::  unorm(ndim),utang(ndim)
	real(prcn) ::  sumforcepoint(3),sumforcenode(3)
	real(prcn) ::  pnmag,dfmag,unmago,unmagb,unmagi, normal(ndim)
	real(prcn) ::tempfor1(ndim),  tempt, tempfor2(ndim),unmagl,ucarl,unorml(ndim),utangl(ndim),ucarf(ndim), tempdist
	real(prcn) ::  drm,utangm(ndim),unormm(ndim),unmag, tmp_rnp, upimod, force_tmp(ndim), ubnorm(ndim), &
		& ubtang(ndim), unormi(ndim), utangi(ndim), ubnmag, utemp
	integer :: isp, i,j,k,l,m,n, is(ndim),iii(ndim),io(ndim), c,d,rcount,bcount,ib, ie, jb, je, kb, ke, onew,&
		& ii, jj, kk , pcell(3), mxftemp
	character*100 :: FILENAME1, FILENAME2, FILENAME3, FILENAME4, FILENAME5, FILENAME6, FILENAME7
	logical:: filexist, isopen, REVERSE
	real(prcn) :: ulb(ndim),plb,plo,pli,deriv_p,deriv_u,deriv_un,frmeanfluidloc(ndim), frmeanloc(ndim), norm_factor, mpg_without_unst(ndim)
	real(prcn), allocatable :: norm_factor_phase(:)

	!real(prcn), ALLOCATABLE,DIMENSION(:,:,:) :: fromleft,fromright
	real(prcn) :: unsteady_term(3)
	real(prcn) :: visc_totalloc(ndim), pres_totalloc(ndim)

	real(prcn) :: l2_norm_du_loc, l2_norm_du
	real(prcn), allocatable :: inner_vel_old(:,:,:)

	Public :: bcset, calc_pgrad, calc_visc, compute_omega, calc_pres_visc_drag
CONTAINS

#if 0
	subroutine compute_mpg_at_n2(rks)
		implicit none 
		integer, Intent(in) ::  rks
		integer ::  n , idim, m, l, iphs
		real(prcn) :: force_total(ndim), x

		if (I_AM_NODE_ZERO)write (*,*)'CALCULATING MPG @ NTH TIME-STEP'
!		call compute_omega
		visc_total_old = visc_total
		pres_total = zero 
		visc_total = zero 
#if PARALLEL
		pres_totalloc = zero 
		visc_totalloc = zero 
#endif

		do m=1,nbody
			if (rks.eq.1) then 
				do n=1,ndim
					pres(m,n)=zero
					visc(m,n)=zero
					torq(m,n)=zero
#if PARALLEL
					presloc(m,n)=zero
					viscloc(m,n)=zero
					torqloc(m,n)=zero
					force_loc(m,n) = zero
#endif
				enddo
			endif
       
			call calc_pres_visc_drag2(m,rks)
		enddo

		do idim=1,ndim
			pres_avg(idim) = SUM(pres(1:nbody,idim))!/real(nspec1,prcn)
			visc_avg(idim) = SUM(visc(1:nbody,idim))!/real(nspec1,prcn)
		enddo
    
		pres_drag= DSQRT(dot_product(pres_avg(1:ndim), pres_avg(1:ndim)))
		visc_drag = DSQRT(dot_product(visc_avg(1:ndim), visc_avg(1:ndim)))

		pres_drag = pres_drag/norm_factor
		visc_drag = visc_drag/norm_factor
		total_drag = (pres_drag+visc_drag)

		if (FROM_POST) goto 1000

		if (.not.set_mpg) then 
			if (move_particles) then
				mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + cf*(usmean(:)-usmean_des(:)) 
				mpg(:) = mpg(:) + (pres_total(:) -visc_total(:)) * (one/(rhof*(one-mean_volfrac)) + one/(mean_volfrac*rhos))/voldom
				mpg(:) = mpg(:)/(one/rhof - one/rhos)
				frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf*(usmean_des(:)-usmean(:))
			else
				mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + (pres_total(:)/voldom -visc_total(:)/voldom)/(one-maxvolfrac)  
			endif

			mpg_without_unst(:) = pres_total(:)/voldom -(visc_total(:))/voldom  

			if (include_frmeanfluid)  then 
				write (*,*)'INCLUDING FRMEAN FLUID'
				mpg(:) = mpg(:) + frmeanfluid(:)
			endif
			mpg_without_unst(:) = mpg_without_unst(:)/(one - maxvolfrac)
			mpg_without_unst(:) = mpg_without_unst(:)/(coef(rks,1) + coef(rks,2))
		else
			if (move_particles) then
				frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf *(usmean_des(:)-usmean(:))
			endif
		endif
1000	continue 
    
		unsteady_term(:) = cf*(ufmean_des(:)-ufmean(:))*(one-maxvolfrac)
		force_total = zero 

		do m = 1, nbody
			force(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:)) -(coef(rks,1)+coef(rks,2))*mpg(:)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0 
			force_total(:) = force_total(:) + force(m,:)

			force_chem(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:))
		enddo

		if (I_AM_NODE_ZERO) then
			if (move_particles) then
				write (*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom/norm_factor
				write (*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)*norm_factor)
		       
				write (*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)*norm_factor) 

				write (*,'(A30,3(2x,g17.8))') 'MPG TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom/norm_factor
			else
				write (*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom/(one-maxvolfrac)/norm_factor

				write (*,'(A30,3(2x,g17.8))') 'PRES TERbcsetM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)/((one-maxvolfrac))/norm_factor

				write (*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)/((one-maxvolfrac))/norm_factor 

				write (*,'(A30,3(2x,g17.8))') 'MPG TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom/norm_factor
				write (*,'(A30,3(2x,g17.8))') 'MPG WITHOUT UNST = ', (coef(rks,1)+coef(rks,2))*mpg_without_unst(:)*voldom/norm_factor

				write (*,'(A30,3(2x,g17.8))') 'FORCE_TOTAL = ', (coef(rks,1)+coef(rks,2))*force_total(:)/norm_factor

			endif
		endif
		return 
	end subroutine compute_mpg_at_n2


	subroutine calc_pres_visc_drag2(m,rks)
		use mypost_process
		implicit none
		integer, intent(in) :: m,rks

		integer :: i, j, k, ii, jj, kk
		integer :: istart, iend, jstart, jend, kstart, kend
		integer :: local_solid, local_fluid, ibody
		real(prcn) :: local_length, local_ufmean(ndim), local_usmean(ndim), local_volfrac, local_Rem, part_force, local_meanslip(ndim), local_meanslip_mag


		!-----------------------------------------------------------------------


		if (lybyd_small<lybyd-small_number) then
			local_length = lybyd_small*(2*radbdy(1))

			istart = xc(m,1) - local_length/2
			iend   = xc(m,1) + local_length/2
			jstart = xc(m,2) - local_length/2
			jend   = xc(m,2) + local_length/2
			kstart = xc(m,3) - local_length/2
			kend   = xc(m,3) + local_length/2

			local_solid = 0
			local_fluid = 0
			local_ufmean = zero
			do kk=kstart, kend
				do jj=jstart, jend
					do ii=istart, iend
						if (ii<1) then
							i=ii+nx
						elseif (ii>nx) then
							i=ii-nx
						else
							i=ii
						endif

						if (jj<1) then
							j=jj+my
						elseif (jj>my) then
							j=jj-my
						else
							j=jj
						endif

						if (kk<1) then
							k=kk+mz
						elseif (kk>mz) then
							k=kk-mz
						else
							k=kk
						endif

						if (fluid_atijk(i,j,k)) then
							local_fluid = local_fluid+1
							local_ufmean(:) = local_ufmean(:) + ur(i,j,k,:)
						else
							local_solid = local_solid+1
							call find_body(i,j,k,ibody, doml(2)/2)
							local_usmean(:) = local_usmean(:) + velbdy(ibody,:)
						endif
					enddo
				enddo
			enddo

			local_volfrac = dble(local_solid) / dble(local_fluid+local_solid)
			local_ufmean(:) = local_ufmean(:) / local_fluid
		else
			local_volfrac = maxvolfrac
			local_ufmean(:) = ufmean(:)
		endif

		local_usmean(:) = velbdy(m,:)

		local_meanslip(:) = local_ufmean(:) - local_usmean(:)
		local_meanslip_mag = sqrt(dot_product(local_meanslip, local_meanslip))

		local_Rem = (1-local_volfrac) * local_meanslip_mag * (radbdy(m)*2*dx) / vis

		call compute_ibm_drag(local_volfrac, local_Rem, part_force)

!write (*,"(1i,4d15.7)") m, local_meanslip(:), local_meanslip_mag
!write (*,"(3d15.7)") local_volfrac, local_Rem, part_force
!write (*,*)

		part_force = part_force * 3 * pi * (radbdy(m)*2*dx) * vis * (1-local_volfrac)
		visc(m,:) = part_force * local_meanslip(:) + (coef(rks,1)+coef(rks,2)) * mpg(:) * pi * (two*radbdy(m)*dx)**3.d0 /6.d0 
		pres(m,:) = zero

		visc_total(:) = visc_total(:) + visc(m,:)
		pres_total(:) = pres_total(:) + pres(m,:)

	end subroutine calc_pres_visc_drag2

#endif









	subroutine compute_mpg_at_n(rks)
		implicit none 

		integer, Intent(in) ::  rks
		integer ::  n , idim, m, l, iphs
		real(prcn) :: force_total(ndim), x

		if (I_AM_NODE_ZERO) write (*,*)'CALCULATING MPG @ NTH TIME-STEP'

		call compute_omega

		visc_total_old = visc_total

		!pres_total = zero 
		!visc_total = zero 

		pres_totalloc = zero 
		visc_totalloc = zero 

		if (rks.eq.1) then 
			!pres(m,:)=zero
			!visc(m,:)=zero
			!torq(m,:)=zero

			do m=1, nbody
				do idim=1, ndim
					presloc(m,idim) = zero
					viscloc(m,idim) = zero
					torqloc(m,idim) = zero
					force_loc(m,idim) = zero
				enddo
			enddo
		endif

		do m=1,nbody            ! loop over bodies
			if (myid_particles(m)) then
				!if (I_AM_NODE_ZERO) write (*,*) "particle = ", m
     
				iphs = 1!part_array(m)%iphs
				nbnd = phase_array(iphs)%nbnd
				NULLifY(bndarray)
				bndarray => phase_array(iphs)%bndpts

				da(1)=4.*pi*(radbdy(m)*dx)**2./real(nbnd,prcn)
				call calc_pres_visc_drag(m,rks)
			endif
		enddo !CLOSE LOOP OVER ALL BODIES


		!do m=1, nbody
			!do n=1, ndim
				GLOBAL_DOUBLE_SUM(presloc, pres, nbody*ndim, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(viscloc, visc, nbody*ndim, comm_cart_2d)
				if (move_particles) GLOBAL_DOUBLE_SUM(torqloc, torq,nbody*ndim, comm_cart_2d)
			!enddo
		!enddo


		GLOBAL_DOUBLE_SUM(pres_totalloc, pres_total, ndim, comm_cart_2d)
		GLOBAL_DOUBLE_SUM(visc_totalloc, visc_total, ndim, comm_cart_2d)
		!do n=1, ndim
		!	GLOBAL_DOUBLE_SUM(pres_totalloc(n),pres_total(n),1,comm_cart_2d)
		!	GLOBAL_DOUBLE_SUM(visc_totalloc(n),visc_total(n),1,comm_cart_2d)
		!enddo

		do idim=1,ndim
			pres_avg(idim) = SUM(pres(1:nbody,idim))!/real(nspec1,prcn)
			visc_avg(idim) = SUM(visc(1:nbody,idim))!/real(nspec1,prcn)
		enddo

		pres_drag = DSQRT(dot_product(pres_avg(1:ndim), pres_avg(1:ndim)))
		visc_drag = DSQRT(dot_product(visc_avg(1:ndim), visc_avg(1:ndim)))

		pres_drag = pres_drag/norm_factor !/(one-maxvolfrac)
		visc_drag = visc_drag/norm_factor !/(one-maxvolfrac)
		total_drag = (pres_drag+visc_drag)
    
		if (.not.set_mpg) then
			if (move_particles) then
				mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + cf*(usmean(:)-usmean_des(:)) 
				!!$mpg(:) = mpg(:) + (pres_total(:) -(visc_total(:)))*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/voldom
				mpg(:) = mpg(:) + (pres_total(:) -(visc_total(:)))*(one/(rhof*(one-mean_volfrac)) + one/(mean_volfrac*rhos))/voldom
				mpg(:) = mpg(:)/(one/rhof - one/rhos)
				frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf *(usmean_des(:)-usmean(:))
			else
				mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + (pres_total(:)/voldom -(visc_total(:))/voldom)/(one-maxvolfrac)  
			endif

			mpg_without_unst(:) = pres_total(:)/voldom -(visc_total(:))/voldom  

			if (include_frmeanfluid) then 
				write (*,*) 'INCLUDING FRMEAN FLUID'
				mpg(:) = mpg(:) + frmeanfluid(:)
			endif
			!if (.not.move_particles) 
			!mpg(:) = mpg(:)/(one - maxvolfrac)
			mpg_without_unst(:) = mpg_without_unst(:)/(one - maxvolfrac)
			! mpg(:) = mpg(:)/(coef(rks,1) + coef(rks,2))
			mpg_without_unst(:) = mpg_without_unst(:)/(coef(rks,1) + coef(rks,2))
		else
			if (move_particles) then
				frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf *(usmean_des(:)-usmean(:))
			endif
		endif
1000	continue

		unsteady_term(:) = cf*(ufmean_des(:)-ufmean(:))*(one-maxvolfrac)
		force_total = zero 
    
		do m = 1, nbody
			force(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:)) -(coef(rks,1)+coef(rks,2))*mpg(:)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0 
			force_total(:) = force_total(:) + force(m,:)
			force_chem(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:))
		enddo
    
		if (I_AM_NODE_ZERO) then
			if (move_particles) then
				write (*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom /norm_factor

				write (*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:) *(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)) *norm_factor

				write (*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:) *(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)) *norm_factor

				write (*,'(A30,3(2x,g17.8))') 'mpg TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom /norm_factor
			else
				write (*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:) *voldom/(one-maxvolfrac) /norm_factor

				write (*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)/((one-maxvolfrac)) /norm_factor

				write (*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)/((one-maxvolfrac)) /norm_factor 

				write (*,'(A30,3(2x,g17.8))') 'mpg TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom /norm_factor
				write (*,'(A30,3(2x,g17.8))') 'mpg WITHOUT UNST = ', (coef(rks,1)+coef(rks,2))*mpg_without_unst(:)*voldom /norm_factor

				write (*,'(A30,3(2x,g17.8))') 'FORCE_TOTAL = ', (coef(rks,1)+coef(rks,2))*force_total(:) /norm_factor
			endif
		endif
		return 
	end subroutine compute_mpg_at_n
  

	subroutine bcset(rks)!(ur,pr,nlr,onlr)
		use general_funcs, only : instant_file_opener
use mypost_process
		implicit none 

		integer, Intent(in) ::  rks
		real(prcn) ::  c_drag_st(nbody), c_drag_tr(nbody), norm_drag_poly_spec(nphases), norm_drag_poly
		real(prcn) ::  avg_force(ndim),avg_force_chem(ndim),avg_force_spec(nphases,ndim),avg_force_chem_spec(nphases,ndim), &
				&    tmp_ferror_array(nerr_steps),LHS_UF(ndim), MOD_LHS_UF, cpu0, cpu1, cpu2, cpu3 !, cpu10, cpu11


		integer :: idim, iphs, pstart, pend
		real(prcn) :: slip_diff
		real(prcn) :: frsumloc(ndim), frsum(ndim)
		COMPLEX(prcn) :: fsumloc,fsum, frtemp
		!ag=7.8d-4
		character(LEN=80) :: formfile
		character*100 filel2norm

		integer :: unitl2norm


		if (I_AM_NODE_ZERO) call CPU_TIME (CPU2) 


		formfile='formatted'

		FILENAME1 = TRIM(RUN_NAME)//'_norm_drag_chem'
		FILENAME4 = TRIM(RUN_NAME)//'_norm_drag'
		FILENAME6 = TRIM(RUN_NAME)//'_drag_components'
		FILENAME7 = TRIM(RUN_NAME)//'_normdrag_poly'

		if (.not.allocated(norm_factor_phase)) allocate(norm_factor_phase(nphases))
		if (zero_slip) then
			norm_factor = (3.d0*pi*vis*(ucharmod+SMALL_NUMBER)*char_length) * (one-maxvolfrac)
			norm_factor_phase(:) = (3.d0*pi*vis*(ucharmod+SMALL_NUMBER)*phase_array(:)%dia) * (one-maxvolfrac)
		else
			norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)
			norm_factor_phase(:) = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*phase_array(:)%dia)
		endif

		if (TRIM(input_type).eq.'lubtest') norm_factor = (3.d0*pi*vis*(one-maxvolfrac)*char_length)*real(nbody,prcn)
		if (.not.only_dem) call calc_visc

		if (rks.eq.1) then
			force(:,:)=zero
			!zero the drag on the body --> the total, viscous and pressure drags.
		endif

		!call compute_mpg_at_n(rks)

!if (I_AM_NODE_ZERO) call CPU_TIME (CPU10) 

		if (.not.only_dem) then
			if (use_drag_law) then
				!call compute_mpg_at_n2(rks)
			else
				call compute_mpg_at_n(rks)
			endif
		endif

!if (I_AM_NODE_ZERO) call CPU_TIME (CPU11)
!if (I_AM_NODE_ZERO) mpg_time = cpu11-cpu10





		if (move_particles) then
			if (I_AM_NODE_ZERO) call calc_part_statistics(rks)

			if (I_AM_NODE_ZERO) call CPU_TIME (CPU0) 

			call des_time_march(.false.)
			if (I_AM_NODE_ZERO) then
				call CPU_TIME (CPU1)
				dem_time_current = cpu1-cpu0
				dem_time = dem_time + dem_time_current
			endif




			call grid_nodes_insphere
			do m = 1, nbody
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
		endif
    
!		if (FROM_POST) then
!			if (use_drag_law) then
!				call compute_mpg_at_n2(rks)
!			else
!				call compute_mpg_at_n(rks)
!			endif
!			call calc_local_pres_visc_drag_plane
!		endif

		if (only_dem) return

		drm = 1.0
		radm = 1.0

		furf = 1.0

!if (I_AM_NODE_ZERO)  call CPU_TIME (CPU10)

		call calc_pgrad

!if (I_AM_NODE_ZERO)  call CPU_TIME (CPU11)
!if (I_AM_NODE_ZERO) pgrad_time = cpu11-cpu10


		!ALLOCATE(fromleft(my,mz,ndim),fromright(my,mz,ndim))
		sumforcepoint(:)=zero
		sumforcenode(:)=zero
    
		apmax = 0.d0
		norm_drag_spec = 0.d0
		norm_drag = 0.d0
		norm_drag_chem_spec = zero
		norm_drag_chem = zero
		norm_drag_poly = zero
		norm_drag_poly_spec = zero

		frmean(:) = zero
		frmeanloc(:) = zero
		frmeanfluid(:) = zero
		frmeanfluidloc(:) = zero
		fr(:,:,:,:)=zero
		!-----------------------------------------------------------------------
		!     loop over all bodies

		l2_norm_du_loc = zero
		l2_norm_du = zero

!if (I_AM_NODE_ZERO) call CPU_TIME (CPU10)

		do m=1,nbody            ! loop over bodies
			if (myid_particles(m)) then
     			!---------------------------------------------------------------------------
        		!     zero the drag on the body --> the total, viscous and pressure drags.
	        	iphs = 1!part_array(m)%iphs
				nbnd = phase_array(iphs)%nbnd
     			nrpr = phase_array(iphs)%nrpr
        		NULLifY(bndarray)
	        	bndarray => phase_array(iphs)%bndpts

				da(1) = 4.*pi*(radbdy(m) *dx)**2./real(nbnd,prcn)
				da(2) = 4.*pi*(radibdy(m)*dx)**2./real(nbnd,prcn) !/real(part_array(m)%nrpr_active, prcn)

	        	!if (dobnd) call calc_bnd_forcing(m,rks)
				if (dorpr) call calc_inner_reversal_forcing(m,rks)
!#if PARALLEL
				!!$if (set_umean) then
				!!$do n = 1, ndim
				!!$GLOBAL_doUBLE_SUM(force_loc(m,n),force(m,n),1,comm_cart_2d)
				!!$enddo
				!!$endif
!#endif
			endif
		enddo

#if 0
		if (debug_check) then
			GLOBAL_DOUBLE_SUM(l2_norm_du_loc,l2_norm_du,1,comm_cart_2d)
			l2_norm_du = sqrt(l2_norm_du)

			if (iglobstep>1) then
				if (I_AM_NODE_ZERO) then
					unitl2norm = 1
					filel2norm = trim(run_name)//"_du_l2norm"
					call instant_file_opener(filel2norm,unitl2norm,.true.)
					write (unitl2norm,'(500(2x, e20.12))') t/t_conv, l2_norm_du
					close(unitl2norm)
				endif
			endif
		endif
#endif

		if (iturbon.and.forced_turbulence) then
			if (t>activate_forcing_time*eddy_time_i) then
				if (I_AM_NODE_ZERO) write (*,*) "LINEAR_FORCING ACTIVE...."
				call linear_forcing
			endif
		endif

    
		GLOBAL_DOUBLE_SUM(frmeanloc,frmean,3,comm_cart_2d)
		if (include_frmeanfluid) GLOBAL_DOUBLE_SUM(frmeanfluidloc,frmeanfluid,3,comm_cart_2d)
		!do idim=1, ndim
		!	GLOBAL_DOUBLE_SUM(frmeanloc(idim),frmean(idim),1,comm_cart_2d)
		!	if (include_frmeanfluid) GLOBAL_DOUBLE_SUM(frmeanfluidloc(idim),frmeanfluid(idim),1,comm_cart_2d)
		!enddo

		!new frmean(1:ndim) = frmean(1:ndim)/real(mx1*my*mz,prcn)
		frmean(1:ndim) = frmean(1:ndim) / real(global_n(1)*global_n(2)*global_n(3),prcn)
		frmeanfluid(1:ndim) = frmeanfluid(1:ndim) / real(count_fluid,prcn)

		call communicate_ib_forces

		if (debug_check) then
			do idim=1, ndim
				frsumloc(idim) = SUM(fr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))
				!GLOBAL_DOUBLE_SUM(frsumloc(idim),frsum(idim),1,comm_cart_2d)
			enddo
			GLOBAL_DOUBLE_SUM(frsumloc,frsum,3,comm_cart_2d)

			if (I_AM_NODE_ZERO) then
				call screen_separator(40,'^')
				write (*,*) "CHECKING FRMEAN:"
				write (*,'(1A,3(2x,g17.8))') 'FRMEAN(1:ndim)=               ', frmean(1:ndim)
				write (*,'(1A,3(2x,g17.8))') 'FRMEAN(1:NDIM) FROM  DOMAIN = ', frsum(1:ndim)/(global_n(1)*global_n(2)*global_n(3))
				call screen_separator(40,'-')
			endif
		endif

		do idim = 1, ndim 
			fr(:,:,:,idim) = fr(:,:,:,idim) - frmean(idim)
			call fftwr2c(fr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim), ff(:,:,:,idim))
		enddo


!if (I_AM_NODE_ZERO)  call CPU_TIME (CPU11)
!if (I_AM_NODE_ZERO) reversal_time = cpu11-cpu10

		!-----------------------------------------------------------------------
		!      POST PROCESSING AND DIAGNOSTICS
		!-----------------------------------------------------------------------
		pstart = 1
		do iphs = 1, nphases
			pend = pstart + phase_array(iphs)%npart - 1
			do idim=1,ndim
				avg_force_spec(iphs,idim) = SUM(force(pstart:pend,idim))/phase_array(iphs)%npart
				avg_force_chem_spec(iphs,idim) = SUM(force_chem(pstart:pend,idim))/phase_array(iphs)%npart
			enddo
			pstart = pend + 1
		enddo

		avg_force = zero
		avg_force_chem = zero

		do iphs = 1, nphases
			avg_force(1:ndim) = avg_force(1:ndim) + phase_array(iphs)%volfrac*avg_force_spec(iphs,1:ndim)
			avg_force_chem(1:ndim) = avg_force_chem(1:ndim) + phase_array(iphs)%volfrac*avg_force_chem_spec(iphs,1:ndim)
		enddo
    
		avg_force(1:ndim) = avg_force(1:ndim)/mean_volfrac
		avg_force_chem(1:ndim) = avg_force_chem(1:ndim)/mean_volfrac

		!!$do iphs = 1, nphases
		!!$if (I_AM_NODE_ZERO)write (*,'(A,4(2x,g17.8))') 'FORCES: avg_force_spec : ', phase_array(iphs)%volfrac,avg_force_spec(iphs,1:ndim)
		!!$enddo
		!!$if (I_AM_NODE_ZERO)write (*,'(A,3(2x,g17.8))') 'FORCES: avg_force : ', avg_force(1:ndim)
		do iphs = 1, nphases
			norm_drag_spec(iphs) = DSQRT(avg_force_spec(iphs,1)**2.d0 + avg_force_spec(iphs,2)**2.d0 + avg_force_spec(iphs,3)**2.d0)
			norm_drag_chem_spec(iphs) = DSQRT(avg_force_chem_spec(iphs,1)**2.d0 + avg_force_chem_spec(iphs,2)**2.d0 + avg_force_chem_spec(iphs,3)**2.d0)
		enddo
    
		norm_drag = DSQRT(avg_force(1)**2.d0 + avg_force(2)**2.d0 + avg_force(3)**2.d0)
		norm_drag_chem = DSQRT(avg_force_chem(1)**2.d0 + avg_force_chem(2)**2.d0 + avg_force_chem(3)**2.d0)
		mpg_avg = DSQRT(dot_product(mpg(1:ndim), mpg(1:ndim)))
		total_drag_mpg = mpg_avg*voldom /norm_factor/nbody
		LHS_UF(:) = dufmeandt(:)*voldom -pres_total(:)/((one-maxvolfrac)) + visc_total(:)/((one-maxvolfrac))
		MOD_LHS_UF = DSQRT(dot_product(LHS_UF(1:3), LHS_UF(1:3)))
		MOD_LHS_UF = MOD_LHS_UF /norm_factor /nbody
		if (I_AM_NODE_ZERO) write (*,'(A,5(2x,g17.8))') 'DRAG: PRES, VISC, TOTAL, MOD_LHS_UF/N, TOTAL_MPG', pres_drag, visc_drag, total_drag, MOD_LHS_UF, total_drag_mpg 
    
		if (TRIM(input_type).eq.'lubtest') then
			do iphs = 1, nphases
				norm_drag_spec(iphs) = norm_drag_spec(iphs) / norm_factor_phase(iphs)
				norm_drag_poly_spec(iphs) = norm_drag_chem_spec(iphs)
				norm_drag_chem_spec(iphs) = norm_drag_chem_spec(iphs) / norm_factor_phase(iphs)
				norm_drag_poly_spec(iphs) = norm_drag_poly_spec(iphs)/(vis**2.d0)
			enddo
			norm_drag  = norm_drag / norm_factor
			norm_drag_poly  = norm_drag_chem
			norm_drag_chem  = norm_drag_chem / norm_factor
			norm_drag_poly  = norm_drag_poly /(vis**2.d0)
		else
			do iphs = 1, nphases
				norm_drag_spec(iphs) = norm_drag_spec(iphs) / norm_factor_phase(iphs)
				norm_drag_poly_spec(iphs) = norm_drag_chem_spec(iphs)
				norm_drag_chem_spec(iphs) = norm_drag_chem_spec(iphs) / norm_factor_phase(iphs)
				norm_drag_poly_spec(iphs) = norm_drag_poly_spec(iphs) / norm_factor_phase(iphs)
			enddo
			norm_drag  = norm_drag / norm_factor
			norm_drag_poly  = norm_drag_chem
			norm_drag_chem  = norm_drag_chem / norm_factor
			norm_drag_poly  = norm_drag_poly / (vis**2.d0)
		endif

		do iphs = 1, nphases
			if (norm_drag_spec(iphs).gt.ZERO) then
				phase_array(iphs)%ferror = ABS(norm_drag_spec(iphs) - phase_array(iphs)%fold) /norm_drag_spec(iphs)
			else
				phase_array(iphs)%ferror = ONE
			endif
		enddo

		if (norm_drag.gt.ZERO) then
			ferror = ABS(norm_drag - fold) / norm_drag
		else
			ferror = one
		endif
		fold = norm_drag
    
		do iphs = 1, nphases
			phase_array(iphs)%fold = norm_drag_spec(iphs)
		enddo

		!!$if (Re.EQ.Zero.and.ReT.gt.zero) then
		!!$norm_drag1 = (norm_drag1)/(6.d0*pi*vis*SQRT(gran_temp)*radbdy(1)*dx)
		!!$else
		!!$norm_drag1 = (norm_drag1)/(6.d0*pi*vis*ucharmod*radbdy(1)*dx)
		!!$norm_drag2 = norm_drag2/(6.d0*pi*vis*ucharmod*radbdy(nbody)*dx)
		!!$norm_drag  = norm_drag/(3.d0*pi*vis*ucharmod*dia_phys)
		!!$endif

		if (rks.eq.itrmax) then 
			!Rearrange the ferror_array array so that the last entry is flushed out
			tmp_ferror_array(1:nerr_steps) = ferror_array(1:nerr_steps)
			ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
			ferror_array(1) = ferror
			!!$PRINT*,'FERROR_A =', FERROR_ARRAY
			ferror_hist = SUM(ferror_array(1:nerr_steps))/nerr_steps
			do iphs = 1, nphases
				tmp_ferror_array(1:nerr_steps) = phase_array(iphs)%ferror_array(1:nerr_steps)
				phase_array(iphs)%ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
				phase_array(iphs)%ferror_array(1) = phase_array(iphs)%ferror
				!!$PRINT*,'FERROR_A =', FERROR_ARRAY
				phase_array(iphs)%ferror_hist = SUM(phase_array(iphs)%ferror_array(1:nerr_steps))/nerr_steps
			enddo
		endif

		if (I_AM_NODE_ZERO) then
			write (*,'(A25,4(2x,g17.8))') 'NORM DRAGS, FERROR', norm_drag, ferror
			if (.not.((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'psd'))) write (*,'(A25,4(2x,g17.8))') &
					& 'NORM DRAGS PHASES:', norm_drag_spec(1:nphases)
			if (rks.eq.itrmax) then 
				!c_drag_st(1:nbody) = force(1:nbody,1)/(3.d0*pi*vis*dia_phys*ucharmod)

				unitnormdragchem = 1
				call instant_file_opener(filename1,unitnormdragchem,.true.)
				write (unitnormdragchem,'(500(2x, e20.12))') t/t_conv, t/t_vis, &
							&	t/t_diff, norm_drag_chem, norm_drag_chem_spec(1:nphases)
				close(unitnormdragchem)	

				unitnormdrag = 1
				call instant_file_opener(FILENAME4,unitnormdrag,.true.)
				write (unitnormdrag,'(500(2x, e20.12))') t/t_conv, t/t_vis, t/t_diff, norm_drag,& 
					&	norm_drag_spec(1:nphases), usmean_act(1:ndim), &
					&	(ufmean(idim)/ucharmod,idim=1,ndim)
				close(unitnormdrag)


				unitdrag_comps = 1
				call instant_file_opener(FILENAME6,unitdrag_comps,.true.)
				write (unitdrag_comps,'(20(2x, e20.12))') pres_drag, visc_drag, total_drag, &
					&	total_drag_mpg, norm_drag, abs(total_drag-total_drag_mpg)/(total_drag + SMALL_NUMBER)
				close(unitdrag_comps)
          
				unitnormdragpoly = 1
				call instant_file_opener(FILENAME7,unitnormdragpoly,.true.)
				write (unitnormdragpoly,'(500(2x, e20.12))') t/t_conv, t/t_vis, t/t_diff, &
					&	norm_drag_poly, norm_drag_poly_spec(1:nphases)
				close(unitnormdragpoly)
			endif
		endif

		! Checking for blown-up simulations
		if (norm_drag.gt.1E+06) then 
			if (I_AM_NODE_ZERO) then 
				write (*,'(A,2x,g17.8,2x,A)') 'NORM DRAG', norm_drag,' is greater than 1E+06'
				write (*,'(A)') 'STOPING THIS CASE AFTER WRITING THE BLOW UP INDICATOR FILE'
				open(2000, file=TRIM(RUN_NAME)//'_BLOWUP.dat', form='formatted')
				close(2000, status="keep")
			endif
			PARALLEL_FINISH()
			stop
		endif

1000	FORMAT(18(E14.6,1x))
2000	continue

		if (I_AM_NODE_ZERO) then
			call CPU_TIME (CPU3) 
			bc_time_current = cpu3-cpu2 - dem_time_current
			bc_time = bc_time + bc_time_current
		endif
		return
	end subroutine bcset

#if 0
  subroutine calc_bnd_forcing(m,rks)
    implicit none
    integer, Intent(in) :: m,rks
    integer :: l

    real(prcn) ::  xl(ndim),xlo(ndim),xli(ndim)
    real(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
    real(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, df(nbnd,ndim), CROSSP(ndim),linvel(ndim)
    logical :: velterm, pressterm
    integer :: count, vcell(ndim)
#if PARALLEL
    real(prcn) :: xltemp
    integer :: bndcount,bndtot, vcelltemp, pcelltemp,focus_point,focus_particle
    
    focus_point = -1
    focus_particle = -1
#endif
    !-----------------------------------------------------------------------
    !write (*,*) 'IN FLOW BND'
    bcount = 0
    frombnd = .true.
    fromrpr = .false.
    count = 0
    
    BNDLOOP: do l=1, nbnd

       rad = zero
       do n=1,ndim

          !     note: x co-ordinate of xc() is absolute center (center with
          !     respect to global origin) minus foffset
          !     Hence, is(1) is the currently the cell location w.r.t foffset

          xl(n)=xc(m,n)+bndarray(n,l)*radbdy(m)

          rad=rad+(bndarray(n,l)*radbdy(m))**2.0


          is(n)=INT(xl(n))
          
          ul(n)=zero
          nll(n)=zero
          onll(n)=zero

          ppll(n)=zero
          dfll(n)=zero
       enddo
       rad = DSQRT(rad)
       xpb(1) = xl(1)-0.5
       xpb(2:3)=xl(2:3)
       do n = 1, ndim
          if (xpb(n).lt.zero) then 
             pcell(n) = int(xpb(n)-1)
             !because of int of example -1.1 is -1, but we want -2. So,
             ! the adjustment 
          else 
             pcell(n) = int(xpb(n))
          endif
          if (xl(n).lt.zero) then
             vcell(n) = int(xl(n)-1)
          else 
             vcell(n) = int(xl(n))
          endif
       enddo

#if PARALLEL
       vcelltemp = vcell(1)
       xltemp  = xl(1)
       if (l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE) then
          PRINT*,' INNER REV PT = ', xli(1),vcelltemp, myid, m
          PRINT*,' EXTERNAL REV PT = ', xlo(1), myid,m
       endif
       
       if (.not.CELL_IN_PROC(vcelltemp)) then
          WEST_PERIODIC_IMAGE(vcell(1),vcelltemp,xl(1),xltemp)
          EAST_PERIODIC_IMAGE_MOD(vcell(1),vcelltemp, xl(1), xltemp)
          if (.not.CELL_IN_PROC(vcelltemp)) then
             if (l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE) then
                PRINT*,' INNER REVERSAL PT = ', xl(1),vcelltemp, myid, m
                
             endif

             goto 2600
          endif
       endif

       velterm = CELL_IN_VEL_GRID(vcelltemp)
       
       vcell(1) = vcelltemp
       xl(1) = xltemp

       pressterm = .true.
       pcelltemp = pcell(1)
       xltemp = xpb(1)
       if (.not.CELL_IN_PRESS_GRID(pcelltemp)) then
          WEST_PERIODIC_IMAGE_PRES(pcell(1),pcelltemp, xpb(1), xltemp)
          EAST_PERIODIC_IMAGE_PRES(pcell(1),pcelltemp, xpb(1), xltemp)
          if (l.eq.FOCUS_POINT) then
             PRINT*,' INNER PRES PT = ', xpi(1),pcelltemp, myid, m
          endif
          pressterm = CELL_IN_PRESS_GRID(pcelltemp)
       endif
       pcell(1) = pcelltemp
       xpb(1) = xltemp
#else
       velterm = .true.
       pressterm = .true.
#endif
       if (pressterm) then

          pl=zero
          call interpolate_pdata(pcell,xpb,ppll,pl,l)

       endif

       call interpolate_udata(vcell,xl,ib&
            &,ie,jb,je,kb,ke,ul,nll,onll,dfll, 1,m, l, onew) 
       
       call CROSS_PRODUCT3(CROSSP(1:3),ANGV(M,1:3,1),SNORM(1:3))
       linvel(1:ndim) = velbdy(m,1:ndim) + CROSSP(1:ndim)*RADBDY(m)*dx 

       
       do n=1,3
          force_tmp(n) = zero
          if (velterm)force_tmp(n) = cf*(-linvel(n)+ ul(n))
          
          force_tmp(n) = force_tmp(n) - coef(rks,3)*nll(n)-coef(rks,4)&
               &*onll(n)
          force_tmp(n) = force_tmp(n) + (coef(rks,1)+coef(rks,2))*(ppll(n)-dfll(n))
          
          force_tmp(n) = force_tmp(n)*da(1)*drm*dx

          sumforcepoint(n)=sumforcepoint(n)+force_tmp(n)
!!$#if PARALLEL          
!!$          if (set_umean) force_loc(m,n)= force_loc(m,n) - force_tmp(n)
!!$#else
!!$          if (set_umean) force(m,n)= force(m,n) - force_tmp(n)
!!$#endif
       enddo

              
       do k = 1, onew 
          do j = 1, onew
             do i = 1, onew

                ii = ib+i-1
                jj = jb+j-1
                kk = kb+k-1
#if !PARALLEL
                if (ii.lt.1) ii = mxf+ii-1
                if (ii.gt.mxf-1) ii = ii-(mxf-1)
#endif
                if (jj.lt.1) jj = my+jj
                if (jj.gt.my) jj = jj-my
                if (kk.lt.1) kk = mz+kk
                if (kk.gt.mz) kk = kk-mz 

                LOCAL_INDEX(ii)

                do n=1,ndim
                   tmppa = weightp(i,j,k)*force_tmp(n)
                   fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0)
                   if (fluid_atijk(ii,jj,kk)) then
                      if (include_frmeanfluid)  frmeanfluid(n) = frmeanfluid(n) + tmppa/(dx**3.d0)
                   endif
#if PARALLEL
                   frmeanloc(n) = frmeanloc(n) + tmppa/(dx**3.d0)
#else
                   frmean(n) = frmean(n) + tmppa/(dx**3.d0)
#endif
                                      
#if 0
						if ((GLOBAL_INDEX(ii).eq.1)) then
							if ((pcell(1).eq.my/2)) then
								fromleft(jj,kk,n) = fromleft(jj,kk,n) + tmppa/(dx**3.d0)
							elseif (pcell(1).eq.my/2+1) then
								fromright(jj,kk,n) = fromright(jj,kk,n) + tmppa/(dx**3.d0)
							endif
						endif
#endif
						enddo
					enddo
				enddo
			enddo
#if PARALLEL
2600      continue
#endif          

		enddo BNDLOOP
		frombnd = .false.
	end subroutine calc_bnd_forcing
#endif

	subroutine calc_inner_reversal_forcing(m,rks)
		implicit none
		integer, Intent(in) :: m,rks
		integer :: l, count_fl, count_so, idim
		!integer :: pcelli(ndim)
		integer :: vcelli(ndim), vcellitmp(ndim)!, vcelli_select(ndim)
		integer :: vcello(ndim), vcellotmp(ndim)
		real(prcn) :: xli(ndim), xlitmp(ndim) !, xli_select(ndim)
		real(prcn) :: xlo(ndim), xlotmp(ndim)
		!real(prcn) :: xpi(ndim) !, force_tmp_vis(ndim)
		real(prcn) :: force_fl(ndim), force_dist(ndim)
		real(prcn) :: ulo(ndim), uli(ndim), tmppa, xltemp,CROSSP(ndim),linvel(ndim)
		!logical :: velterm, pressterm
		logical :: i_have_the_inner_cell, i_have_the_outer_point, find_boundary_cell

		!write (*,*) 'IN FLOW REVERSAL, nrpr :', nrpr
		frombnd = .false.
		fromrpr = .true.


#if 0
		if (debug_check) then
			if (.not.allocated(inner_vel_old)) then
				allocate(inner_vel_old(nbody,nrpr,ndim))
				inner_vel_old = zero
			endif
		endif
#endif


		do l=1,nrpr
			if (.not.PART_ARRAY(M)%if_rev(L)) goto 666

			!     location of internal points
			xli(:)=xc(m,:)+ bndarray(:,l)*radibdy(m)

#if 0
if (abs(bndarray(1,l))<small_number .and. abs(bndarray(2,l))<small_number .and. bndarray(3,l)>small_number ) then
	check_point=.true.
else
	check_point=.false.
endif

if (check_point) then
	write (*,"(1a,3d15.7)") "bndarray = ", bndarray(:,l)
	write (*,"(1a,3d15.7)") "xc       = ", xc(l,:)
	read (*,*)
endif
#endif

			uli(:)=zero
			onll(:)=zero
			nll(:)=zero
			ppll(:)=zero
			dfll(:)=zero

			!     location of external points
			xlo(:)=xc(m,:)+ bndarray(:,l)*radobdy(m)
			ulo(:)=zero

			rad= sqrt(dot_product(bndarray(:,l),bndarray(:,l))) * radibdy(m)

			plb = zero
			plo = zero
			pli = zero
			snorm(:) = (bndarray(:,l)*radibdy(m))/rad

			do n = 1, ndim
				vcelli(n) = floor(xli(n))
				vcello(n) = floor(xlo(n))
			enddo

			i_have_the_inner_cell = cell_in_this_domain(vcelli,vcellitmp,xli,xlitmp)
			if (i_have_the_inner_cell) then
				i_have_the_outer_point = point_in_this_ghost_domain(vcello, vcellotmp, xlo, xlotmp)

				!THIS LOOP IS TO MAKE SURE THE INNER AND OUTER POINTS THAT ARE MUTUAL
				!AMONG ADJACENT DOMAINS ARE ONLY CHOSEN BY ONE OF THE DOMAINS
				find_boundary_cell = .false.
				do idim=1, ndim
					if (vcellitmp(idim)==0 .or. vcellitmp(idim)==local_ni(idim)) then
						!vcelli_select(:) = vcellitmp(:)
						!xli_select(:)    = xlitmp(:)
						find_boundary_cell = .true.
						exit
					endif
				enddo

				if (find_boundary_cell .and. i_have_the_outer_point) then
					do idim=1, ndim
						!THIS "IF" SHOULD BE EXTENDED ALSO TO "IDIM==1" IF THE DOMAIN IS
						!ALSO DECOMPOSED ALONG THE X-DIRECTION
						if ((idim==2.and.nprocy>1) .or. (idim==3.and.nprocz>1)) then
							if (vcellitmp(idim)==0) then
								if (xli(idim) <= xlo(idim)) then
								else
									i_have_the_outer_point = .false.
									exit
								endif
							elseif (vcellitmp(idim)==local_ni(idim)) then
								if (xli(idim) > xlo(idim)) then
								else
									i_have_the_outer_point = .false.
									exit
								endif
							endif
						endif
					enddo
				endif
			else
				i_have_the_outer_point = .false.
			endif		

			if (debug_check) then
				if (in_how_many_domains(i_have_the_outer_point)/=1) then
					if (I_AM_NODE_ZERO) then
						write (*,"(1a,8i6)") "THIS OUTER POINT IN NO DOMAIN, OR MORE THAN ONE DOMAIN, M, l : ", m, l, vcello, vcellotmp
					endif
					if (i_have_the_outer_point) write (*,"(1a,7i4,6d15.7)") "DOMAIN, INDEX, POSITION ", myid, vcellitmp, vcellotmp, xlitmp, xlotmp

					PARALLEL_FINISH()
					stop
				endif
			endif

			if (i_have_the_outer_point) then
				vcelli(:) = vcellitmp(:)
				xli(:)    = xlitmp(:)

				vcello(:) = vcellotmp(:)
				xlo(:)    = xlotmp(:)

				!-----------------------------------------------------------------------
				!     calculate outer local velocities
				!-------------------------------------------------------
				call interpolate_udata(vcello,xlo,ib,ie,jb,je,kb,ke,ulo,nll,onll,dfll, 0, m, l, onew) 


				!unmago=zero
				!do n=1,ndim
				!	unmago = unmago + snorm(n)*ulo(n)
				!enddo

				!-----------------------------------------------------------------------
				!     calculate internal local pressure and velocities
				!-------------------------------------------------------
				call interpolate_pdata(vcelli,xli,ppll,pl,l)
				call interpolate_udata(vcelli,xli,ib,ie,jb,je,kb,ke,uli,nll,onll,dfll, 1, m, l, onew) 

				!unmagi = dot_product(snorm(:),ulo(:))

				unmag  = zero
				ubnmag = zero! Velocity of the body (dot) normal vector
				call CROSS_PRODUCT3(CROSSP(1:3),ANGV(M,1:3,1),SNORM(1:3))
				linvel(1:ndim) = velbdy(m,1:ndim) + CROSSP(1:ndim)*RADBDY(m)*dx 

				unmag  = dot_product(snorm(:), ulo(:))
				ubnmag = dot_product(snorm(:), linvel(:))

				!-----------------------------------------------------------------------
				!     set internal forcing to reverse external velocity
				!     reverse both tangential velocity components
				!     zero radial component
				!     scale tangential velocities by ratio of radii
				!-----------------------------------------------------------------------
				ucar(:)=zero
				do d=1,ndim
					unorm(d)=snorm(d)*unmag
					utang(d)=ulo(d)-unorm(d)
					ubnorm(d) = snorm(d)*ubnmag  ! velocity of the body in the normal direction
					!ubtang(d) = velbdy(m,d)-ubnorm(d) ! velocity of the body in the tangential direction
					ubtang(d) = linvel(d)-ubnorm(d) ! velocity of the body in the tangential direction
					utangi(d) = ubtang(d)*(radobdy(m)-radibdy(m))-utang(d)*(radbdy(m)-radibdy(m))
					utangi(d) = utangi(d)/(radobdy(m)-radbdy(m))

					! Find tangential velocity at the internal reversal point so that the no slip condition is satisfied
					if (revernorm.EQ.1) then
						unormi(d) = ubnorm(d)*(radobdy(m)-radibdy(m))-unorm(d)*(radbdy(m)-radibdy(m))
						unormi(d) = unormi(d)/(radobdy(m)-radbdy(m))
						!unormi(d) = 2*ubnorm(d) * (radbdy(m)/radibdy(m))**2 - unorm(d) * (radobdy(m)/radibdy(m))**2
					else
						unormi(d) = zero
					endif

					if (.not.bubble_particles) then
						ucar(d) =  utangi(d) + unormi(d)
					else
						ucar(d) = -utangi(d) + unormi(d)
					endif
				enddo

				do n=1,ndim
					force_tmp(n) = zero
					force_tmp(n)=cf*(-ucar(n)+uli(n))
	       
					force_tmp(n) = force_tmp(n) - coef(rks,3)*nll(n)-coef(rks,4)*onll(n)
					force_tmp(n) = force_tmp(n) + (coef(rks,1)+coef(rks,2))*(ppll(n)-dfll(n))
					force_tmp(n) = force_tmp(n)*da(2)*drm*dx
					!force_tmp_vis(n) = force_tmp_vis(n)*da(1)*drm*dx
					sumforcepoint(n) = sumforcepoint(n)+force_tmp(n)
				enddo

#if 0
				if (debug_check) then
					if (iglobstep>1) then
						l2_norm_du_loc = l2_norm_du_loc + dot_product(uli(:)-inner_vel_old(m,l,:),uli(:)-inner_vel_old(m,l,:))
					endif
					inner_vel_old(m,l,:) = ucar(:)
				endif
#endif

				count_fl = 0
				force_fl = zero
				do k = 1, onew
					kk = kb+k-1
					do j = 1, onew
						jj = jb+j-1
						do i = 1, onew
							ii = ib+i-1
							if (fluid_atijk(ii,jj,kk)) then 
								count_fl = count_fl+1
								do n=1,ndim
									tmppa = weightp(i,j,k)*force_tmp(n)
									force_fl(n) = force_fl(n) + tmppa/(dx**3.d0)
								enddo
							endif
						enddo
					enddo
				enddo

				count_so = onew*onew*onew - count_fl
				force_dist(:) = force_fl(:)/real(count_so)

				do k = 1, onew
					kk = kb+k-1
					do j = 1, onew
						jj = jb+j-1
						do i = 1, onew
							ii = ib+i-1
							do n=1,ndim                   
								tmppa = weightp(i,j,k)*force_tmp(n)
								if (include_frmeanfluid)  then 
									if (fluid_atijk(ii,jj,kk)) then
										frmeanfluidloc(n) = frmeanfluidloc(n) + tmppa/(dx**3.d0)
									endif
									fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0)
								else
									if (.not.fluid_atijk(ii,jj,kk)) then 
										fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0) + force_dist(n)
									endif
								endif

								frmeanloc(n) = frmeanloc(n) + tmppa / (dx**3.d0)
							enddo
						enddo
					enddo
				enddo
			endif


#if 0
if (check_point) then
	write (*,"(1a,3i)") "vcello = ", vcello(:)
	write (*,"(1a,3i)") "vcelli = ", vcelli(:)

	write (*,"(1a,3d15.7)") "ulo   = ", ulo(:)
	write (*,"(1a,3d15.7)") "uli   = ", uli(:)

	write (*,"(1a,3d15.7)") "snorm = ", snorm(:)
	write (*,"(1a,3d15.7)") "ubody = ", linvel

	write (*,"(1a,3d15.7)") "unamg = ", unmag
	write (*,"(1a,3d15.7)") "unorm = ", unorm
	write (*,"(1a,3d15.7)") "utang = ", utang

	write (*,"(1a,3d15.7)") "unormi = ", unormi
	write (*,"(1a,3d15.7)") "utangi = ", utangi
	write (*,"(1a,3d15.7)") "ucar   = ", ucar

	write (*,"(1a,3d15.7)") "fr_tmp = ", force_tmp
	write (*,"(1a,3d15.7)") "fr_tmp = ", force_fl
	read (*,*)
endif
#endif



			!------------------------------------------
			!     Close loop over all reversal points
666		continue
88		enddo		! loop over reversal points

		fromrpr = .false.
	end subroutine calc_inner_reversal_forcing

	subroutine calc_pres_visc_drag(m,rks)
		USE bcsetarrays, ONLY :  omega => fr
		implicit none
		integer, Intent(in) :: m,rks
		integer :: l, pcelltemp,vcelltemp, pcellb(ndim), vcellb(ndim)

		real(prcn) ::  xl(ndim),xlo(ndim),xli(ndim)
		real(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
		real(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, df(nbnd,ndim)
		real(prcn) :: tempforce(ndim), crossp(ndim), xltemp, xptemp

		logical :: i_have_the_point
		integer :: index_out(ndim)
		real(prcn) :: position_out(ndim)
		!-----------------------------------------------------------------------
		!write (*,*) 'IN FLOW BND'
		bcount = 0
		da(1)=4.*pi*(radbdy(m)*dx)**2./real(nbnd,prcn)

		do l=1,nbnd
			if (.not.part_array(m)%if_drag(L)) goto 100
			rad = zero
			do n=1,ndim
				xl(n)=xc(m,n)+ bndarray(n,l)*radbdy(m)

				rad=rad+(bndarray(n,l)*radbdy(m))**2.0
				is(n)=INT(xl(n))

				ul(n)=zero
				ppll(n)=zero
       	enddo

			rad = DSQRT(rad)

			xpb(1) = xl(1) !new -0.5
			xpb(2:3)=xl(2:3)
			do n = 1, ndim
				vcellb(n) = floor(xl(n))
				pcellb(n) = floor(xpb(n))

				!if (xpb(n).lt.zero) then 
				!	pcellb(n) = int(xpb(n)-1)
				!	!because of int of example -1.1 is -1, but we want -2. So,
				!	! the adjustment 
				!else 
				!	pcellb(n) = int(xpb(n))
				!endif
				!if (xl(n).lt.zero) then 
				!	vcellb(n) = int(xl(n)-1)
				!else 
				!	vcellb(n) = int(xl(n))
				!endif
			enddo

			i_have_the_point = point_in_this_domain(vcellb, index_out, xl, position_out)

			if (debug_check) then
				if (in_how_many_domains(i_have_the_point)/=1) then
					if (I_AM_NODE_ZERO) write (*,"(1a,6i10)") "THIS BOUNDARY NODE IN NO DOMAIN, OR MORE THAN ONE DOMAIN : ", vcellb, index_out
					if (i_have_the_point) write (*,"(1a,7i4,6d15.7)") "DOMAIN, INDEX, POSITION ", myid, vcellb, index_out, xl, position_out

					PARALLEL_FINISH()
					stop
				endif
			endif

			if (i_have_the_point) then
				vcellb(:) = index_out(:)
				xl(:) = position_out(:)

				pcellb(:) = vcellb(:)
				xpb(:) = xl(:)

				pl=zero
				call interpolate_pdata(pcellb,xpb,ppll,pl,l)
				call interpolate_udata(vcellb,xl,ib,ie,jb,je,kb,ke,ul,nll,onll,dfll, 0, m, l, onew) 

				vort(:) = zero 
				do k = 1, onew 
					do j = 1, onew
						do i = 1, onew
							ii = ib+i-1
							jj = jb+j-1
							kk = kb+k-1
!#if !PARALLEL
							!if (ii.lt.1) ii = mx+ii !new ii = mxf+ii-1
							!if (ii.gt.mx) ii = ii-mx !new if (ii.gt.mxf-1) ii = ii-(mxf-1)
!#endif
							!if (jj.lt.1) jj = my+jj
							!if (jj.gt.my) jj = jj-my
							!if (kk.lt.1) kk = mz+kk
							!if (kk.gt.mz) kk = kk-mz 
							!LOCAL_INDEX(ii)                

							do n=1, ndim
								vort(n)=vort(n)+ weightp(i,j,k)*omega(ii,jj,kk,n)
							enddo
						enddo
					enddo
				enddo
				df(l,1)=vort(3)*cd(2,l)
				df(l,1)=df(l,1)-vort(2)*cd(3,l)

				df(l,2)=vort(1)*cd(3,l)
				df(l,2)=df(l,2)-vort(3)*cd(1,l)

				df(l,3)=vort(2)*cd(1,l)
				df(l,3)=df(l,3)-vort(1)*cd(2,l)

				df(l,1)=-df(l,1)
				df(l,2)=-df(l,2)
				df(l,3)=-df(l,3)
				!---------------------------------------------------------------
				!     calculate the viscous and pressure components separately
				!---------------------------------------------------------------
				do d=1, ndim
					presloc(m,d)= presloc(m,d)-pl*cd(d,l)*da(1)
					pres_totalloc(d) = pres_totalloc(d) + pl*cd(d,l)*da(1)

					viscloc(m,d)=viscloc(m,d)+df(l,d)*vis*da(1)
					visc_totalloc(d) = visc_totalloc(d) + df(l,d)*vis*da(1)

					tempforce(d) = df(l,d)*vis*da(1)!+ df(l,d)*vis*da(1)
				enddo
				call CROSS_PRODUCT3(CROSSP(1:3),bndarray(1:3,l),TEMPFORCE(1:3))

				TORQLOC(m,1:3) = TORQLOC(m,1:3) + CROSSP(1:3)*RADBDY(m)*dx
				!---------------------------------------------
				!     close loop over all boundary points
				!---------------------------------------------
			endif
100		continue
		enddo
	end subroutine calc_pres_visc_drag


	subroutine compute_omega
		USE bcsetarrays, ONLY :  omega => fr
		implicit none
		integer i, j, k, dim1, dim2, dim3
		complex(prcn) :: wtmp
		real(prcn) :: epsilon

		do dim1=1, ndim
			uftmp(:,:,:) = czero
			do dim2=1, ndim
				if (dim2==dim1) goto 10
				do dim3=1, ndim
					if (dim3==dim1.or.dim3==dim2) goto 20
#if !PARALLEL
					do k=1, local_no(3)
						do j=1, local_no(2)
							do i=1, local_no(1)
								if (dim2==1) then
									wtmp = wx(i)
								elseif (dim2==2) then
									wtmp = wy(j)
								elseif (dim2==3) then
									wtmp = wz(k)
								endif
#else
					do k=1, local_no(2)
						do j=1, local_no(1)
							do i=1, local_no(3)
								if (dim2==1) then
									wtmp = wy(j)
								elseif (dim2==2) then
									wtmp = wz(k)
								elseif (dim2==3) then
									wtmp = wx(i)
								endif
#endif
								if (dim1==1 .and. dim2==2 .and. dim3==3) then
									epsilon = one
								elseif (dim1==3 .and. dim2==1 .and. dim3==2) then
									epsilon = one
								elseif (dim1==2 .and. dim2==3 .and. dim3==1) then
									epsilon = one
								elseif (dim1==dim2 .or. dim1==dim3 .or. dim2==dim3) then
									epsilon = zero
								else
									epsilon = -one
								endif

								uftmp(i,j,k) = uftmp(i,j,k) + epsilon * u(i,j,k,dim3) * wtmp
							enddo
						enddo
					enddo
20					continue
				enddo
10				continue
			enddo
			! 3d c_to_r fftw uftmp -> omega(:,:,:,dim1)
			call fftwc2r(uftmp, omega(1:local_ni(1),1:local_ni(2),1:local_ni(3),dim1))
			call communicate_in_gohst_domain(omega(:,:,:,dim1))
		enddo
	end subroutine compute_omega

	subroutine calc_visc
		implicit none
		integer :: idim, i, j, k

		do idim=1, ndim
#if !PARALLEL
			do k=1, local_no(3)
				do j=1, local_no(2)
					do i=1, local_no(1)
#else
				do k=1, local_no(2)
					do j=1, local_no(1)
						do i=1, local_no(3)
#endif
						uftmp(i,j,k) = -w2(i,j,k)*u(i,j,k,idim)*vis
					enddo
				enddo
			enddo
			! 3d c_to_r fftw uftmp -> diffn
			call fftwc2r(uftmp, diffn(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))

			call communicate_in_gohst_domain(diffn(:,:,:,idim))
		enddo
	end subroutine calc_visc
     
	subroutine calc_pgrad
		implicit none
		integer :: i,j,k,idim
		complex(prcn) :: wtmp
       
       !c     To Check the contribution of pressure terms in the forcing
       !c     Compute pressure gradient and store at --pressure-- grid points
       !c     earlier was being stored at velocity grid points, and implicit
       !c     smoothing of the pressure field was being done in the x-direction

		do idim=1, ndim
#if !PARALLEL
			do k=1, local_no(3) !mz
				do j=1, local_no(2) !my
					do i=1, local_no(1) !mx2
						if (idim==1) then
							wtmp = wx(i)
						elseif (idim==2) then
							wtmp = wy(j)
						elseif (idim==3) then
							wtmp = wz(k)
						endif
#else
			do k=1, local_no(2)
				do j=1, local_no(1)
					do i=1, local_no(3)
						if (idim==1) then
							wtmp = wy(j)
						elseif (idim==2) then
							wtmp = wz(k)
						elseif (idim==3) then
							wtmp = wx(i)
						endif
#endif
						uftmp(i,j,k) = wtmp*p(i,j,k)
					enddo
				enddo
			enddo

			call fftwc2r(uftmp, ppr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))

			ppr(:,:,:,idim) = ppr(:,:,:,idim) + mpg(idim)

			call communicate_in_gohst_domain(ppr(:,:,:,idim))
		enddo
	end subroutine calc_pgrad

#if 0     
     subroutine calc_local_pres_visc_drag_plane
       USE bcsetarrays, ONLY :  omega => fr
       implicit none
       logical, SAVE :: routine_called = .false. 
       integer :: m
       integer :: l, iphi, itheta, unitno, unitno1, unitno2, unitno3, unitno4, unitno5
       real(prcn) ::  xl(ndim),xlo(ndim),xli(ndim)
       real(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
       real(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, df(nbnd,ndim)
       real(prcn) :: tempforce(ndim), crossp(ndim), pres_loc, visc_loc,&
            & pres_bnd(3), visc_bnd(3), force_bdy_tmp(3), force_bdy_mag,&
            & axial_direction(3), thetaang, theta
       real(prcn) :: cphi, phiang, xcor, xsmin, xsmax, xmin, xmax, dx_x,&
            & nx, ny,nz, rad_proj, dtheta, normal(3), total_Fvis(nbody&
            &,2), total_Pres(nbody,2), areabdy(nbody), norm_factor2,&
            & avg_Fvis(2), avg_Pres(2), Fvis_theta_avg(count_theta),&
            & Fvis_phi_avg(count_phi),  Pres_phi_avg(count_phi),&
            & Pres_theta_avg(count_theta), conf1, conf2
       integer :: norm_phi(count_phi), norm_theta(count_theta)

       norm_factor2 = (3.d0*pi*vis*(meanslipmod+SMALL_NUMBER)*dia_phys) 
       total_Fvis = zero 
       total_pres = zero 

       write (*,*) 'count_phi = ', count_phi, count_theta
       axial_direction(1:3) = uchar(1:3)/ucharmod
       unitno = getnewunit(minunitno, maxunitno)

       open(unitno,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_PHI_THETAZERO.dat', form="form&
            &atted",status="unknown") 
       unitno1 = getnewunit(minunitno, maxunitno)
       open(unitno1,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_THETA_PHIZERO.dat', form="form&
            &atted",status="unknown") 

       unitno2 = getnewunit(minunitno, maxunitno)
       open(unitno2,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_PHI_NBODY.dat', form="form&
            &atted",status="unknown") 

       unitno3 = getnewunit(minunitno, maxunitno)
       open(unitno3,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_THETA_NBODY.dat', form="fo&
            &rmatted",status="unknown")  


       unitno4 = getnewunit(minunitno, maxunitno)
       open(unitno4,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_PHI_AVG.dat', form="fo&
            &rmatted",status="unknown") 

       unitno5 = getnewunit(minunitno, maxunitno)
       open(unitno5,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_THETA_AVG.dat', form="fo&
            &rmatted",status="unknown")  

       dtheta  = twopi/real(count_theta-1,prcn)

       NULLifY(bndarray)
       bndarray => phase_array(1)%bndpts

       xsmin = MINVAL(bndarray(1,1:nbnd))
       xsmax = MAXVAL(bndarray(1,1:nbnd))

       if (.not.routine_called) then 
          ALLOCATE(Fvis_theta(nbody, count_theta), Fvis_phi(nbody,&
               & count_phi), Pres_theta(nbody, count_theta), Pres_phi(nbody, count_phi))
          routine_called = .true.
       endif
       L = 1
       force_bdy_tmp = zero 
       !-----------------------------------------------------------------------
       write (*,*) 'CALCULATING FVIS AND PRES DRAG COMPONENTS ALONG THETA &
            &AND PHI' 

       BODYLOOP: do m = 1, nbody !Loop over bodies  
          areabdy(m) = 4.*pi*(radbdy(m)*dx)**2
          norm_phi = zero 
          norm_theta = zero 
          Fvis_phi(m,:) = zero 
          Fvis_theta(m,:) = zero 
          Pres_phi(m,:) = zero 
          Pres_theta(m,:) = zero 

          write (unitno,*)'Zone'
          write (unitno1,*)'Zone'
          write (unitno2,*)'Zone'
          write (unitno3,*)'Zone'

          xmin = xsmin!*radbdy(m)
          xmax = xsmax!*radbdy(m)

          dx_x = (xmax-xmin)/real(count_phi-one,prcn)

          PHILOOP: do iphi = 1, count_phi
             nx = xmin + dx_x*real((iphi-1),prcn)

             xl(1)=xc(m,1)+nx*radbdy(m) !global coordinate

             cphi=(nx)

             phiang = ACOS(cphi)
             rad_proj = sin(phiang)
             phiang = oneeighty*phiang/pi
             !write (*,'(A,5(2x,g17.8))')'xcor= ', iphi, nx, phiang, rad_proj
             phi_array(iphi) = phiang
             THETALOOP: do itheta = 1, count_theta 
                theta = dtheta*real(itheta-one,prcn)
                thetaang = oneeighty*theta/pi
                !if (iphi.eq.1) write (*,*) 'theta =', thetaang
                ny = rad_proj*cos(theta)
                nz = rad_proj*sin(theta) 
                
                xl(2)=xc(m,2)+ny*radbdy(m) !global coordinate
                xl(3)=xc(m,3)+nz*radbdy(m) !global coordinate
                
                theta_array(itheta) = theta*oneeighty/pi
                norm_phi(iphi) = norm_phi(iphi) + 1
                norm_theta(itheta) = norm_theta(itheta) + 1
                
                !print*,iphi, itheta, norm_phi(iphi), norm_theta(itheta)
                do n=1,ndim
                   
                   is(n)=INT(xl(n))
                   
                   ul(n)=zero

                   ppll(n)=zero

                enddo
                rad = radbdy(m)

                pl=zero
                isp=INT(xl(1)-0.5)
                xpb(1) = xl(1)-0.5
                xpb(2:3)=xl(2:3)
                do n = 1, ndim
                   if (xpb(n).lt.zero) then 
                      pcell(n) = int(xpb(n)-1)
                      !because of int of example -1.1 is -1, but we want -2. So,
                      ! the adjustment 
                   else 
                      pcell(n) = int(xpb(n))
                   endif
                enddo


                call interpolate_pdata(pcell,xpb, ppll,pl,l)

                do n = 1, ndim
                   if (xl(n).lt.zero) then 
                      pcell(n) = int(xl(n)-1)
                   else 
                      pcell(n) = int(xl(n))
                   endif
                enddo

                call interpolate_udata(pcell,xl,ib&
                     &,ie,jb,je,kb,ke,ul,nll,onll,dfll, 0,m, l, onew) 

                vort(:) = zero 
                do k = 1, onew 
                   do j = 1, onew
                      do i = 1, onew
                         ii = ib+i-1
                         jj = jb+j-1
                         kk = kb+k-1
                         if (ii.lt.1) ii = mxf+ii-1
                         if (ii.gt.mxf-1) ii = ii-(mxf-1)
                         if (jj.lt.1) jj = my+jj
                         if (jj.gt.my) jj = jj-my
                         if (kk.lt.1) kk = mz+kk
                         if (kk.gt.mz) kk = kk-mz 

                         do n=1,ndim
                            vort(n)=vort(n)+ weightp(i,j,k)*omega(ii,jj,kk,n) 
                         enddo
                      enddo
                   enddo
                enddo

                normal(1) = nx
                normal(2) = ny
                normal(3) = nz

                !write (*,'(6(2x,g17.8))') 'norm = ', normal(1:3), vort(1:3)
                df(L,1)=vort(3)*normaL(2)
                df(L,1)=df(L,1)-vort(2)*normal(3)

                df(L,2)=vort(1)*normal(3)
                df(L,2)=df(L,2)-vort(3)*normal(1)

                df(L,3)=vort(2)*normal(1)
                df(L,3)=df(L,3)-vort(1)*normal(2)

                df(L,1)=-df(L,1)
                df(L,2)=-df(L,2)
                df(L,3)=-df(L,3)

                !if (THETAANG.GT.oneeighty) write (*,*)'p1 = ', pl
                !---------------------------------------------------------------
                !     calculate the viscous and pressure components separately
                !---------------------------------------------------------------

                do d=1,ndim,1
                   pres_bnd(d) = -pl
                   visc_bnd(d) = df(l,d)
                enddo

                visc_loc = array_dot_product(visc_bnd(1:3),&
                     & axial_direction(1:3))
                pres_loc = -pl*array_dot_product(normal(1:3),&
                     & axial_direction(1:3))

                visc_loc = vis*visc_loc*areabdy(m)/norm_factor2
                pres_loc = pres_loc*areabdy(m)/norm_factor2

                !write (*,'(6(2x,g17.8))') 'norm = ', pres_loc, visc_loc
                !pres_loc = pres_loc/(0.5d0*meanslipmod*meanslipmod)
                !visc_loc = visc_loc*dchar/meanslipmod
                Fvis_phi(m,iphi) = Fvis_phi(m,iphi) + visc_loc

                Fvis_theta(m,itheta) = Fvis_theta(m,itheta) + visc_loc

                pres_phi(m,iphi) = pres_phi(m,iphi) + pres_loc
                pres_theta(m,itheta) = pres_theta(m,itheta) + pres_loc
                
                if (thetaang.eq.zero) then 
                   write (unitno,21) oneeighty-phiang,pres_loc,visc_loc

                endif
                

                if (phiang.eq.half*oneeighty) then 
                   write (unitno1,21) thetaang,pres_loc,visc_loc
                endif
             enddo THETALOOP
             

          enddo PHILOOP
          !READ(*,*)
          do j = 1, count_phi 
             Fvis_phi(m,j) = Fvis_phi(m,j)/real((norm_phi(j)),prcn)
             Pres_phi(m,j) = Pres_phi(m,j)/real((norm_phi(j)),prcn)

             total_Fvis(m,1) =  total_Fvis(m,1) + Fvis_phi(m,j)
             total_pres(m,1) =  total_Pres(m,1) + Pres_phi(m,j)

             write (unitno2,21) oneeighty-phi_array(j), Fvis_phi(m,j),&
                  & pres_phi(m,j) 

          enddo
          total_Fvis(m,1) = total_Fvis(m,1)/count_phi
          total_Pres(m,1) = total_Pres(m,1)/count_phi

          do j = 1, count_theta
             Fvis_theta(m,j) = Fvis_theta(m,j)/real((norm_theta(j)),prcn)
             Pres_theta(m,j) = Pres_theta(m,j)/real((norm_theta(j)),prcn)

             total_Fvis(m,2) = total_Fvis(m,2) + Fvis_theta(m,j)
             total_pres(m,2) = total_Pres(m,2) + Pres_theta(m,j)

             write (unitno3,21) theta_array(j), Fvis_theta(m,j),&
                  & pres_theta(m,j)
          enddo

          total_Fvis(m,2) = total_Fvis(m,2)/count_theta
          total_Pres(m,2) = total_Pres(m,2)/count_theta

          write (*,'(A40,/,2x, i4, 2(2x,g17.8))') 'AVG VISC DRAG  ALONG PHI AND THET&
               &A FOR M = ',M, total_Fvis(m,1), & 
               & total_Fvis(m,2)


          write (*,'(A40,/,2x, i4, 2(2x,g17.8))') 'AVG PRES DRAG  ALONG PHI AND THET&
               &A FOR M = ',M, total_Pres(m,1), & 
               & total_Pres(m,2)

       enddo BODYLOOP

       avg_Fvis(1) =  SUM(total_Fvis(:,1))/real(nbody,prcn)
       avg_Fvis(2) =  SUM(total_Fvis(:,2))/real(nbody,prcn)
       avg_Pres(1) =  SUM(total_Pres(:,1))/real(nbody,prcn)
       avg_Pres(2) =  SUM(total_Pres(:,2))/real(nbody,prcn)

       do j = 1, count_phi 
          Fvis_phi_avg(j) = SUM(Fvis_phi(1:nbody, j))/real(nbody,prcn)
          Pres_phi_avg(j) = SUM(Pres_phi(1:nbody, j))/real(nbody,prcn)
       enddo

       do j=1, count_phi
          conf1 = zero; conf2 = zero
          if (nbody.ge.2) then  
             do m = 1, nbody
                conf1 = conf1 + (Fvis_phi(m,j) - Fvis_phi_avg(j))&
                     &**two
                conf2 = conf2 + (Pres_phi(m,j) - Pres_phi_avg(j))&
                     &**two
             enddo
             conf1 = dsqrt(conf1/real(nbody-1,prcn))
             conf1=1.96/sqrt(dble(nbody))*conf1 
             conf2 = dsqrt(conf2/real(nbody-1,prcn))
             conf2=1.96/sqrt(dble(nbody))*conf1 
          endif
          write (unitno4,21) oneeighty-phi_array(j), Fvis_phi_avg(j),&
               & pres_phi_avg(j) , conf1, conf2 
       enddo

       do j = 1, count_theta
          Fvis_theta_avg(j) = SUM(Fvis_theta(1:nbody, j))/real(nbody,prcn)
          Pres_theta_avg(j) = SUM(Pres_theta(1:nbody, j))/real(nbody,prcn)
       enddo

       do j=1, count_theta
          conf1 = zero; conf2 = zero
          if (nbody.ge.2) then 
             do m = 1, nbody
                conf1 = conf1 + (Fvis_theta(m,j) - Fvis_theta_avg(j))&
                     &**two
                conf2 = conf2 + (Pres_theta(m,j) - Pres_theta_avg(j))&
                     &**two

             enddo
             conf1 = dsqrt(conf1/real(nbody-1,prcn))
             conf1=1.96/sqrt(dble(nbody))*conf1 
             conf2 = dsqrt(conf2/real(nbody-1,prcn))
             conf2=1.96/sqrt(dble(nbody))*conf1 
          endif

          write (unitno5,21) theta_array(j), Fvis_theta_avg(j),&
               & pres_theta_avg(j) , conf1, conf2 
       enddo

       write (*,'(A40,/, 4(2x,g17.8))') 'AVG VISC DRAG  ALONG PHI AND THET&
            &A  = ', avg_Fvis(1:2), SUM(Fvis_phi_avg(1:count_phi))&
            &/real(count_phi),  SUM(Fvis_theta_avg(1:count_theta))&
            &/real(count_theta)


       write (*,'(A40,/, 4(2x,g17.8))') 'AVG PRES DRAG  ALONG PHI AND THET&
            &A', avg_Pres(1:2), SUM(Pres_phi_avg(1:count_phi))&
            &/real(count_phi),  SUM(Pres_theta_avg(1:count_theta))&
            &/real(count_theta) 
21     FORMAT(10(2xe17.4))

       CLOSE(unitno,status= 'keep')
       CLOSE(unitno1,status= 'keep')
       CLOSE(unitno2,status= 'keep')
       CLOSE(unitno3,status= 'keep')
       CLOSE(unitno4,status= 'keep')
       CLOSE(unitno5,status= 'keep')
     end subroutine calc_local_pres_visc_drag_plane


  subroutine write_complex_forcing
    implicit none
    integer :: i, j, k

    if (I_AM_NODE_ZERO) then
       
       open(1000,FILE=TRIM(RUN_NAME)//'_complex_force.dat',status='unknown')
       write (1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" ', ' "FY" &
            &',' "FZ" '!,' "P" ' !, ' "nl2" ', ' "nl3" '
#if PARALLEL
       write (1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
#else
       write (1000,*)'ZONE F=POINT, I=', nx/2,  ', J=', my, ', K=', mz
#endif
       do k=1,mz
          do j=1,my2
#if PARALLEL
             do i=1,nx !mx
#else
             do  i = 1,nx/2
#endif
                write (1000,*)(GLOBAL_INDEX(i)),(j),(k),dreal(ff(i,j,k,1))&
                           &,dreal(ff(i,j,k,2)),dreal(ff(i,j,k,3))
!!$
!!$                      write (1000,*)(GLOBAL_INDEX(i)),(j),(k),fr(i,j,k,1)&
!!$                           &,fr(i,j,k,2),fr(i,j,k,3)!,pr(i,j,k)/(half*upi(1)*upi(1))!,divur(i,j,k)!,nlbcp(i,j,k,1),nlbcp(i,j,k,2),nlbcp(i,j,k,3)!!velr_mxf(i,j,k,3)
!!$                      !write (14,*)(i),(j),(k),divur(i,j,k)
             enddo
          enddo
       enddo
       close(1000,status='keep')
    endif

    end subroutine write_complex_forcing

  subroutine write_real_forcing
    USE nlmainarrays , ONLY : nlbc, ubc
    implicit none
    integer :: i, j, k
        

#if PARALLEL       
    open(1000,FILE='_real_force.dat',status='unknown')
#else
    open(1000,FILE=TRIM(RUN_NAME)//'_real_force.dat',status='unknown')
#endif
    write (1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" ', ' "FY" &
               &',' "FZ" '
#if PARALLEL
    write (1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
#else
    write (1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
#endif
       do k=1,mz
          do j=1,my
#if PARALLEL
             do i=1,nx !mx
#else             
             do  i =  1, nx
#endif
!             write (1000,*)(GLOBAL_INDEX(i)),(j),(k),(fromleft(j,k,1))&
!                  &,(fromright(j,k,1)),(fromleft(j,k,1)+fromright(j,k,1)), nlbc(i,j,k,1)
                
!             write (1000,*)(GLOBAL_INDEX(i)),(j),(k),(fromleft(j,k,1))&
!                  &,(fromright(j,k,1)),(fromleft(j,k,1)+fromright(j,k,1)), nlbc(i,j,k,1)

             write (1000,*)(GLOBAL_INDEX(i)),(j),(k),(fr(i,j,k,1))&
                  &,(ppr(i,j,k,1)),(nlbc(i,j,k,1))


!!$
!!$                      write (1000,*)(GLOBAL_INDEX(i)),(j),(k),fr(i,j,k,1)&
!!$                           &,fr(i,j,k,2),fr(i,j,k,3)!,pr(i,j,k)/(half*upi(1)*upi(1))!,divur(i,j,k)!,nlbcp(i,j,k,1),nlbcp(i,j,k,2),nlbcp(i,j,k,3)!!velr_mxf(i,j,k,3)
!!$                      !write (14,*)(i),(j),(k),divur(i,j,k)
             enddo
          enddo
       enddo
       close(1000,status='keep')
!       endif
     end subroutine write_real_forcing

#endif

end MODULE boundary_condition



#if 0


			if (i_have_the_inner_cell) then
				i_have_the_outer_cell = cell_in_this_ghost_domain(vcello,vcellotmp,xlo,xlotmp)

#if PARALLEL
				!THIS IS TO MAKE SURE THE INNER AND OUTER POINTS THAT ARE MUTUAL
				!BETWEEN TWO ADJACENT DOMAINS ARE ONLY CHOSEN BY ONE OF THE DOMAINS

				if (i_have_the_outer_cell) then
					do idim=1, ndim
						if (vcelli(idim)==0) then
							index_in(:)    = vcello(:)
							position_in(:) = xlo(:)
							i_have_the_point = point_in_this_domain(vcellb, index_out, xl, position_out)
							


							if (xlitmp(idim)<=xlotmp(idim)) then
							else
								i_have_the_outer_cell = .false.
								exit
							endif
						elseif (vcellitmp(idim)==local_ni(idim)) then
							if (xlotmp(idim)<xlitmp(idim)) then
							else
								i_have_the_outer_cell = .false.
								exit
							endif
						endif
					enddo
				endif




#if 0
				if (i_have_the_outer_cell) then
					do idim=1, ndim
						if (vcellitmp(idim)==0) then
							if (xlitmp(idim)<=xlotmp(idim)) then
							else
								i_have_the_outer_cell = .false.
								exit
							endif
						elseif (vcellitmp(idim)==local_ni(idim)) then
							if (xlotmp(idim)<xlitmp(idim)) then
							else
								i_have_the_outer_cell = .false.
								exit
							endif
						endif
					enddo
				endif
#endif
#endif
			else
				i_have_the_outer_cell = .false.
			endif

			!MAKE SURE THAT THIS FORCING POINT HAS BEEN SELECTED BY ONLY ONE PROCESSOR
			if (debug_check) then
				if (in_how_many_domains(i_have_the_outer_cell)/=1) then
					if (I_AM_NODE_ZERO) then
						write (*,"(2(1A,1I))") "ERROR FOR NODE ", L, " ON PARTICEL ", M
						write (*,"(1a54,3i10,3d15.7)") "AN OUTER NODE IN NO DOMAIN, OR MORE THAN ONE DOMAIN : ", vcello, xlo
						write (*,"(1a54,3i10,3d15.7)") " ", vcellotmp, xlotmp
						write (*,"(1a54,3i10,3d15.7)") "CORRESPONDING TO THE INNEER POINT : ", vcelli, xli
						write (*,"(1a54,3i10,3d15.7)") " ", vcellitmp, xlitmp
					endif
					if (i_have_the_outer_cell) write (*,"(1a,1i,1a)") "DOMAIN ", myid, " HAS THIS DOMAIN."

					PARALLEL_FINISH()
					stop
				endif
			endif
#endif



