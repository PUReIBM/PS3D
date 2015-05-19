module fes

#include "ibm.h"
	use precision 
	use constants 
	use global_data
	use dependent_functions
	use nlmainarrays, ONLY : ur=>ubcp, pr=>pbcp, urstar=>onlbc
	use bcsetarrays, ONLY : ppr
	use nl_allflow
	use usteptime
	use boundary_condition, only : calc_pres_visc_drag, compute_omega
	use report
	use machine
	use surface_module

	integer :: rks, itercount
	real(prcn) :: umean_star(ndim), betak(ndim), afk(ndim)
	real(prcn) :: forcemod, forcemod_old
	logical :: iterstop


contains

	subroutine fes_solver(itr)
		implicit none

		integer :: i, j, k, m, idim
		integer, intent(in) :: itr
		real(prcn) :: cpu0, cpu1
		character*100 surface_file

		rks = itr
		hydroforce_time_current = zero
		fourier_time_current = zero

		if (I_AM_NODE_ZERO) call cpu_time(cpu0) 
		if (iglobstep==1) call calc_velreal(u, umean, ur)		
		call calc_pressure
		call compute_mean_fluid_velocity
		if (I_AM_NODE_ZERO) then
			call cpu_time(cpu1) 
			fourier_time_current = fourier_time_current + cpu1-cpu0
		endif

		call fes_compute_new_timestep(.true.)


		if (move_particles) then
			Write (*,'(A)') 'Computing hydrodynamic forces on particles...'
			if (I_AM_NODE_ZERO) call cpu_time(cpu0) 
			call compute_hydrodynamic_forces
			if (I_AM_NODE_ZERO) then
				call cpu_time(cpu1) 
				hydroforce_time_current = hydroforce_time_current + cpu1-cpu0
			endif
		endif


		Write (*,'(A)') 'Computing convective term using conservative formulation...'
		call fes_nonlinear

		t = t + dt !New time level		

		write (*,'(A)') 'Predicting Fluctuating Velocity...'
		!Predicts velocity with Adams-Bashforth
		! for convective terms, Crank-Nicolson for diffusive terms,
		! pressure at previous time level and no IB forcing term
		if (I_AM_NODE_ZERO) call cpu_time(cpu0) 

		force_factor = zero
		call fes_predict_velocity 
		if (I_AM_NODE_ZERO) then
			call cpu_time(cpu1) 
			vp_time_current = cpu1-cpu0
		endif

		if (I_AM_NODE_ZERO) call cpu_time(cpu0)
		call calc_velreal(onl, umean_star, urstar)
		if (I_AM_NODE_ZERO) then
			call cpu_time(cpu1) 
			fourier_time_current = fourier_time_current + cpu1-cpu0
		endif

		Write(*,'(A)')'Starting inner iterations...'
		do k = 1, mz
			do j = 1, my
				do i = 1, nx
					if(j<=my2) pstar(i,j,k) = p(i,j,k)
					ur(i,j,k,:) = urstar(i,j,k,:) ! Initialize the real velocity with the predicted velocity for iterations
					ppr(i,j,k,:) = zero
				end do
				! I SHOULD CHECK THIS FOR PARALLEL VERSION
				ur(mx,j,k,:) = ur(1,j,k,:)
			end do
		end do

		bc_time_current = zero

		betak(:) = zero
		afk(:) = zero
		forcemod_old = zero
		iterstop = .false.
		itercount = 0
		do while(.not.iterstop)
			itercount = itercount + 1

			if (itercount==1) then
				surface_file = "01"
			elseif (itercount==2) then
				surface_file = "02"
			elseif (itercount==3) then
				surface_file = "03"
			elseif (itercount==4) then
				surface_file = "04"
			elseif (itercount==5) then
				surface_file = "05"
			elseif (itercount==6) then
				surface_file = "06"
			elseif (itercount==7) then
				surface_file = "07"
			elseif (itercount==8) then
				surface_file = "08"
			elseif (itercount==9) then
				surface_file = "09"
			elseif (itercount==10) then
				surface_file = "10"
			elseif (itercount==11) then
				surface_file = "11"
			endif

			surface_file = "FES"//trim(surface_file)
			call surface_field(surface_file)


			if (move_particles) then
				call grid_nodes_insphere

				do m=1, nbody
					call update_nrpr_array(m)
!					call particle_in_proc(m)
				end do
			endif

			call compute_fes_ibforcing

			if (I_AM_NODE_ZERO) call cpu_time(cpu0)
			do k = 1, mz
				do j = 1, my2
					call fes_poisson(j,k)
				end do
			end do
			if (I_AM_NODE_ZERO) then
				call cpu_time(cpu1) 
				vp_time_current = vp_time_current + cpu1-cpu0
			endif

			!call check_fes_divergence(.true.)
			!Write(*,'(A25,2(2x,g12.5))') 'MAX DIVERGENCE (k) = ', divmax


			if (I_AM_NODE_ZERO) call cpu_time(cpu0)
			do idim = 1, ndim
				do i = 1, nx
					call ff2cr(ff(i,:,:,idim),ppr(i,:,:,idim))
				end do
				!ppr(mx,j,k,idim) = ppr(1,j,k,idim)
			end do

			call correct_mean_quantities
!			call calc_velreal(u, umean, ur)
			call calc_pressure
			if (I_AM_NODE_ZERO) then
				call cpu_time(cpu1) 
				fourier_time_current = fourier_time_current + cpu1-cpu0
			endif

			!Write(*,'(A)')'Computing hydrodynamic forces on particles...'

			if (move_particles) then
				if (I_AM_NODE_ZERO) call cpu_time(cpu0)			
				call compute_hydrodynamic_forces
				if (I_AM_NODE_ZERO) then
					call cpu_time(cpu1) 
					hydroforce_time_current = hydroforce_time_current + cpu1-cpu0
				endif
			endif

			if (move_particles.and.nbody>1) then
				!  Write(*,'(A)')'Computing contact forces on particles...'
!				call compute_contact_forces
			end if

			call check_iter_convergence(iterstop)
		end do
    
		mpg(:) = mpg(:) + betak(:)
		frame_accln(:) = frame_accln(:) + afk(:)

		WRITE(*,'(A25,3(2x,g17.8))')'MPG (n) = ', (MPG(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'AF (n) = ', (FRAME_ACCLN(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'<f>(n) = ', (FRMEAN(idim), idim = 1, ndim)

		if(move_particles.and.I_AM_NODE_ZERO) CALL calc_part_statistics(rks)

		if (.not.move_particles) then
			if (I_AM_NODE_ZERO) call cpu_time(cpu0)			
			call compute_hydrodynamic_forces
			if (I_AM_NODE_ZERO) then
				call cpu_time(cpu1) 
				hydroforce_time_current = hydroforce_time_current + cpu1-cpu0
			endif
		endif

		call fes_compute_new_timestep(.true.)
		call report_force


		surface_file = "_final"

		surface_file = "FES"//trim(surface_file)
		call surface_field(surface_file)

write (*,*) "ENTER"
read (*,*)


	end subroutine fes_solver	

	subroutine compute_hydrodynamic_forces
		implicit none

		integer :: m, n, idim, iphs
		real(prcn) :: norm_factor


		norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)*real(nbody,prcn)

		visc_total_old = visc_total
		pres_total = zero 
		visc_total = zero 
#if PARALLEL
		pres_totalloc = zero 
		visc_totalloc = zero 
#endif

		force(:,:)=zero
		pres(:,:)=zero
		visc(:,:)=zero
		torq(:,:)=zero
#if PARALLEL
		presloc(:,:)=zero
		viscloc(:,:)=zero
		torqloc(:,:)=zero
		force_loc(:,:)=zero
#endif

		call compute_omega

		! loop over bodies
		do m=1, nbody
			iphs = 1!part_array(m)%iphs
			nbnd = phase_array(iphs)%nbnd
			NULLIFY(bndarray)
			bndarray => phase_array(iphs)%bndpts

			call calc_pres_visc_drag(m,rks)
#if PARALLEL
			do n = 1, ndim
				GLOBAL_doUBLE_SUM(presloc(m,n),pres(m,n),1,decomp_group)
				GLOBAL_doUBLE_SUM(viscloc(m,n),visc(m,n),1,decomp_group)
				if(move_particles) GLOBAL_doUBLE_SUM(torqloc(m,n),torq(m,n),1,decomp_group)
			enddo
#endif
			force(m,:) = pres(m,:) + visc(m,:)
		enddo !CLOSE LOOP OVER ALL BODIES


#if PARALLEL    
		GLOBAL_doUBLE_SUM(pres_totalloc(1),pres_total(1),3,decomp_group)
		GLOBAL_doUBLE_SUM(visc_totalloc(1),visc_total(1),3,decomp_group)
#endif
    
		do idim=1,ndim
			pres_avg(idim) = SUM(pres(1:nbody,idim))!/real(nspec1,prcn)
			visc_avg(idim) = SUM(visc(1:nbody,idim))!/real(nspec1,prcn)
		end do

		pres_drag= DSQRT(dot_product(pres_avg(1:ndim), pres_avg(1:ndim)))
		visc_drag = DSQRT(dot_product(visc_avg(1:ndim), visc_avg(1:ndim)))
    
		pres_drag = pres_drag/norm_factor
		visc_drag = visc_drag/norm_factor
		total_drag = (pres_drag+visc_drag)
	end subroutine compute_hydrodynamic_forces


	subroutine fes_nonlinear
		use nlmainarrays, Only :  ubc, nlbc, onlbc
		implicit none 

		real(prcn) ::  eold, tvis, nlmean(ndim), nlmean_loc(ndim) 
		complex(prcn) :: nltemp1
		integer i,j,k,l,n,idim
		real(prcn) cpu0, cpu1

		if (I_AM_NODE_ZERO) call cpu_time(cpu0)

		do n=1,ndim
			do i = 1, nx+1
				do j = 1, my2
					do k = 1, mz
						onl(i,j,k,n) = nl(i,j,k,n)
						nl(i,j,k,n) = czero
					end do
				end do
			end do
		end do
    
		if (I_AM_NODE_ZERO) write(*,'(A50)')'CONSERVATIVE: PERFOMING DEL.(UU)...'
		call form_nl

		do n=1,ndim

       !Jamals convention of absorbing the negative
       ! sign into the nl term 
       
#if !PARALLEL          
			nl(1:nx+1,:,:,n) = -nl(1:nx+1,:,:,n)
#else
			nl(1:nx+1,:,:,n) = -nl(1:nx+1,:,:,n)
#endif
		enddo

		if(debug_check)then
			do n=1,ndim
#if PARALLEL       
				nlbc(0,:,:,n) = zero
				do i=1,nx
#else
				do i = 1, nx+1   
#endif
					call ff2cr(nl(i,:,:,n), nlbc(i,:,:,n))
				enddo
#if PARALLEL       
				nlbc(nx+1,:,:,n) = zero
#endif
			enddo
			nlmean(1:ndim) = zero

			do n = 1, ndim 
				nlmean_loc(n) = SUM(nlbc(1:nx,:,:,n))
				GLOBAL_DOUBLE_SUM(nlmean_loc(n),nlmean(n),1,decomp_group)
			end do
			nlmean(:) = nlmean(:)/(mx1*my*mz) 
			if(I_AM_NODE_ZERO) WRITE(*,'(A25,3(2x,g12.5))')'NLMEAN = ', nlmean
		end if

		if (I_AM_NODE_ZERO) then
			call cpu_time(cpu1) 
			nl_time_current = cpu1-cpu0 - fourier_time_nl

			fourier_time_current = fourier_time_current + fourier_time_nl
		endif
	end subroutine fes_nonlinear

	subroutine fes_predict_velocity
		implicit none
		integer :: j, k, rks, idim

		rks = 1
		do k = 1, mz
			do j = 1, my2
				call ustep(rks,j,k)

				do idim = 1, ndim
					onl(1:nx+1,j,k,idim) = u(1:nx+1,j,k,idim) !Storing the predicted velocity u* in onl
				end do
			end do
		end do

		!Predict mean velocity
		umean_star(:) = umean(:) - dt*mpg(:)/rhof - dt*frame_accln(:)
	end subroutine fes_predict_velocity



	SUBROUTINE compute_fes_ibforcing
		USE bcsetarrays, ONLY :  fr, ppr
		IMPLICIT NONE

		INTEGER :: i, j, k, idim, m, iphs
		LOGICAL :: partproc
		real(prcn) :: frmeanloc(ndim), cpu0, cpu1
    
 
		if (I_AM_NODE_ZERO) call cpu_time(cpu0) 

		fr(:,:,:,:) = zero
		frmean(:) = zero
		frmeanloc(:) = zero

		DO m=1,nbody          
			iphs = 1!part_array(m)%iphs
			nbnd = phase_array(iphs)%nbnd
			nrpr = phase_array(iphs)%nrpr
			NULLIFY(bndarray)
			bndarray => phase_array(iphs)%bndpts
#if PARALLEL
			partproc = partinproc(m)
#else 
			partproc = .TRUE.
#endif
			if (partproc) CALL fes_inner_reversal_forcing(m)
		end DO!CLOSE LOOP OVER ALL BODIES
    
#if PARALLEL    
		GLOBAL_DOUBLE_SUM(frmeanloc(1),frmean(1),ndim,decomp_group)
#endif

		frmean(1:ndim) = frmean(1:ndim)/real(count_total,prcn)

		if (I_AM_NODE_ZERO) then
			call cpu_time(cpu1) 
			bc_time_current = bc_time_current + cpu1-cpu0
		endif

		! CHECK THIS FOR PARALLEL
!		CALL communicate_forcing

		if (I_AM_NODE_ZERO) call cpu_time(cpu0) 
		DO idim = 1, ndim 
#if PARALLEL
			! I DON'T KNOW WHY. I SHOULD CHECK IT FOR THE PARALLEL VERSION
			!DO i = 1,nx
			DO i = 2,nx-1
#else
			DO i = 1, nx
#endif
				do j = 1, my 
					do k = 1, mz
						fr(i,j,k,idim) = fr(i,j,k,idim) - frmean(idim)
					end do
				end do
          
				CALL ff2rc(fr(i,1:my,1:mz,idim), ff(i,1:my2,1:mz,idim))
			END DO
		END DO
#if !PARALLEL
		ff(nx+1, :, :,1:ndim) = ff(1,:,:,1:ndim)
#endif
		if (I_AM_NODE_ZERO) then
			call cpu_time(cpu1) 
			fourier_time_current = fourier_time_current + cpu1-cpu0
		endif
	END SUBROUTINE compute_fes_ibforcing



	subroutine fes_inner_reversal_forcing(m)
		use bcsetarrays, ONLY :  fr
		implicit none
		integer, intent(in) :: m
		integer ::  count_fl, count_so
		integer :: i,j,k,l,n, iii(ndim),io(ndim),ib, ie, jb, je, kb, ke, onew,&
		 & ii, jj, kk, d
		integer :: vcelli(ndim), vcello(ndim)
		real(prcn) ::  xl(ndim),xlo(ndim),xli(ndim), force_fl(ndim), force_dist(ndim)
		real(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, xltemp,CROSSP(ndim),linvel(ndim)
		real(prcn) :: pplo(ndim), ppli(ndim), rad, da, snorm(ndim), ucar(ndim), unorm(ndim)
		real(prcn) :: unmago, unmagi, unmag, ubnmag, drm, utangi(ndim)
		real(prcn) :: utang(ndim), ubnorm(ndim), ubtang(ndim), ubtangi(ndim), unormi(ndim), force_tmp(ndim)
		logical :: velterm
#if PARALLEL
		integer :: rprcountloc, rprcount,rprcom,rprcomloc, FOCUS_POint, VELGL, vcelltemp, FOCUS_PARTICLE
#endif

		frombnd = .false.
		fromrpr = .true.
#if PARALLEL
		focus_point = -1
		focus_particle = -1
#endif

		da = 4.*pi*(radibdy(m)*dx)**2./real(part_array(m)%nrpr_active, prcn)
		drm = one
		do l=1,nrpr
			if(.not.PART_ARRAY(M)%if_rev(L)) goto 666

			rad = zero
			do n=1,ndim
				!     location of internal points
				xli(n)=xc(m,n)+ bndarray(n,l)*radibdy(m)

				iii(n)=int(xli(n))
				uli(n)=zero
				ppli(n)=zero
          
				!     location of external points
				xlo(n)=xc(m,n)+ bndarray(n,l)*radobdy(m)

				io(n)=int(xlo(n))
				ulo(n)=zero
				pplo(n)=zero

				rad=rad+(bndarray(n,l)*radibdy(m))**2.0
			enddo
			rad=dsqrt(rad)

			do n=1,ndim                       
				snorm(n)=(bndarray(n,l)*radibdy(m))/rad
			enddo

			do n = 1, ndim
				if(xlo(n).lt.zero) then 
					vcello(n) = int(xlo(n)-1)
				else 
					vcello(n) = int(xlo(n))
				end if

				if(xli(n).lt.zero) then 
					vcelli(n) = int(xli(n)-1)
				else 
					vcelli(n) = int(xli(n))
				end if
			end do

#if PARALLEL
			vcelltemp = vcelli(1)
			xltemp  = xli(1)
			if (l.eq.FOCUS_POint.and.m.eq.FOCUS_PARTICLE) then
				PRint*,' INNER REV PT = ', xli(1),vcelltemp, myid, m
				PRint*,' EXTERNAL REV PT = ', xlo(1), myid,m
			end if
      
			if(.not.CELL_IN_PROC(vcelltemp))then
				WEST_PERIODIC_IMAGE(vcelli(1),vcelltemp,xli(1),xltemp)
				EAST_PERIODIC_IMAGE_MOD(vcelli(1),vcelltemp, xli(1), xltemp)
				if(.not.CELL_IN_PROC(vcelltemp))then
					if(l.eq.FOCUS_POint.and.m.eq.FOCUS_PARTICLE)then
						PRint*,' INNER REVERSAL PT = ', xli(1),vcelltemp, myid, m
						!PARALLEL_FINISH()
						!STOP
					end if
					!             PRint*,' XLTEMP = ', myid, xltemp, xli(1), l
					goto 666
				end if
			end if

			if (EAST_NO_MANS_LAND(vcelli(1)).or.EAST_NO_MANS_LAND(vcelltemp)) then 
				velterm = .not.CONCAVE(xli,1,m)
			elseif (WEST_NO_MANS_LAND(vcelli(1)).or.WEST_NO_MANS_LAND(vcelltemp)) then
				velterm = CONCAVE(xli,1,m)
			else
				velterm = .true.
			end if
			vcelli(1) = vcelltemp
			xli(1) = xltemp
      
			if(velterm)then
				vcelltemp = vcello(1)
				xltemp = xlo(1)
				if(.not.RPR_CELL_IN_PROC(vcelltemp))then
					WEST_PERIODIC_IMAGE(vcello(1),vcelltemp, xlo(1),xltemp)
					EAST_PERIODIC_IMAGE_MOD(vcello(1),vcelltemp,xlo(1),xltemp)
					if(.not.RPR_CELL_IN_PROC(vcelltemp))then
						if(xstart.eq.1)then
							if(vcelltemp.eq.mxf-3)then
								vcelltemp = vcelltemp-(mxf-1)+1
								xltemp = xltemp-(mxf-1)
							endif
						else
							PRint*,' ERROR WITH EXTERNAL POint IN THIS PROCESSOR : ', myid, m, l, xlo(1), vcelltemp,vcello(1),xli(1)
						endif
					endif
				endif
				vcello(1) = vcelltemp
				xlo(1) = xltemp
			endif
      
			if (l.eq.FOCUS_POint.and.m.eq.FOCUS_PARTICLE) then
				print*,' VELTERM = ', l,velterm, myid, m
				PARALLEL_FINISH()
				stop
			endif
#else
			velterm = .true.
#endif
      
			if (velterm) then
				call fes_interpolate_data(vcello, xlo, ib, ie, jb, je, kb, ke, ulo, pplo, .true., m, l, onew) 
				!ulo(:) External velocity at the latest iteration

				unmago=zero
				do n=1,ndim
					unmago = unmago + snorm(n)*ulo(n) 
				enddo
			end if

			call fes_interpolate_data(vcelli,xli,ib,ie,jb,je,kb,ke,uli,ppli, .false., m, l, onew) 
			!FALSE implies interpolating predicted velocity field at the inner reversal points

			if(velterm)then
				unmagi=zero
				do n=1,ndim
					unmagi = unmagi + snorm(n)*ulo(n)
				enddo

				unmag=zero
				ubnmag = zero

				! CHECK THIS FOR FREELY EVOLVING PARTICLES
				! Velocity of the body (dot) normal vector
!				call CROSS_PRODUCT3(CROSSP(1:3),ANGV(M,1:3),SNORM(1:3))
				CROSSP(1:3) = zero

				linvel(1:ndim) = velbdy(m,1:ndim) + CROSSP(1:ndim)*RADBDY(m)*dx 
				do d=1,ndim
					ucar(d)=zero
					unmag=unmag+snorm(d)*ulo(d)
					ubnmag = ubnmag + snorm(d)*linvel(d)
				enddo

				!-----------------------------------------------------------------------
				!     set internal forcing to reverse external velocity
				!     reverse both tangential velocity components
				!     zero radial component
				!     scale tangential velocities by ratio of radii
				!-----------------------------------------------------------------------

				do d=1,ndim
					unorm(d)=snorm(d)*unmag
					utang(d)=ulo(d)-unorm(d)

					ubnorm(d) = snorm(d)*ubnmag  ! velocity of the body in the normal direction
					ubtang(d) = linvel(d)-ubnorm(d) ! velocity of the body in the tangential direction

					utangi(d) = ubtang(d)*(radobdy(m)-radibdy(m))-utang(d)*(radbdy(m)-radibdy(m))
					utangi(d) = utangi(d)/(radobdy(m)-radbdy(m))

					! Find tangential velocity at the internal reversal point so that the no slip condition is satisfied

					IF(revernorm.EQ.1)THEN
						unormi(d) = ubnorm(d)*(radobdy(m)-radibdy(m))-unorm(d)*(radbdy(m)-radibdy(m))
						unormi(d) = unormi(d)/(radobdy(m)-radbdy(m))
					ELSE
						unormi(d) = zero
					endIF
					ucar(d) = utangi(d) + unormi(d)
				enddo
			end if

			do n=1,ndim
				force_tmp(n) = zero

				if(velterm)then
					force_tmp(n) = cf*(-ucar(n)+uli(n)) + (ppli(n)+betak(n))/rhof + afk(n)
					force_tmp(n) = force_tmp(n)*da*drm*dx
				end if
			enddo
      
			count_fl = 0
			force_fl = zero
      
			do k = 1, onew
				kk = kb+k-1
				if(kk.lt.1) kk = mz+kk
				if(kk.gt.mz) kk = kk-mz 
         
				do j = 1, onew
					jj = jb+j-1
					if(jj.lt.1) jj = my+jj
					if(jj.gt.my) jj = jj-my
            
					do i = 1, onew
						ii = ib+i-1
#if !PARALLEL
						if (ii.lt.1) ii = mxf+ii-1
						if (ii.gt.mxf-1) ii = ii-mxf +1
#endif
						LOCAL_INDEX(ii)
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
			!if(count_fl.gt.0)then
				!Write(*,*) ' count_fl gt 0....node = ', myid, count_fl, l, ib,force_fl(1)
			!end if
			force_dist(:) = force_fl(:)/real(count_so)
       
			do k = 1, onew
				kk = kb+k-1
				if (kk.lt.1) kk = mz+kk
				if (kk.gt.mz) kk = kk-mz 

				do j = 1, onew
					jj = jb+j-1
					if (jj.lt.1) jj = my+jj
					if (jj.gt.my) jj = jj-my

					do i = 1, onew
						ii = ib+i-1
#if !PARALLEL
						if(ii.lt.1) ii = mxf+ii-1
						if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
						LOCAL_INDEX(ii)
						do n=1,ndim                   
							tmppa = weightp(i,j,k)*force_tmp(n)
							if(.not.fluid_atijk(ii,jj,kk)) then 
								fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0) + force_dist(n)
							endif
#if PARALLEL
							frmeanloc(n) = frmeanloc(n) + tmppa/(dx**3.d0)
#else
							frmean(n) = frmean(n) + tmppa/(dx**3.d0)
#endif
						enddo
					enddo
				enddo
			enddo
			!     Close loop over all reversal points
			!------------------------------------------
666		CONTINUE             
88		enddo                ! loop over reversal points
		fromrpr = .false.
	end subroutine fes_inner_reversal_forcing
 
	subroutine fes_interpolate_data(pc, pos, ib, ie, jb, je, kb, ke, ul, ppll, uext, ibody,l,onew)
		use general_funcs
		use interpolation
		use bcsetarrays, ONLY :  ppr
		use nlmainarrays, Only : ur=>ubcp, urstar=>onlbc
		implicit none 
		integer, intent(in) :: pc(3)
		integer, intent(in) :: ibody,l
		logical, intent(in) :: uext
		real(prcn), Dimension(:), intent(in) :: pos
		integer, intent(out) :: ib, ie, jb,je,kb,ke, onew
		real(prcn), intent(out), Dimension(:) :: ul, ppll
		integer :: i, j,k, ii,jj,kk, n


		call set_interpolation_stencil(pc, ib, ie, jb, je, kb, ke, interp_scheme, onew) 
		if ((ib.lt.1.or.ie.gt.mxf).and..not.xperiodic) write (*,*) 'Error in i ....',ib,ie,pc,pos
		do k = 1, onew
			do j = 1, onew
				do i = 1, onew
					ii = ib+i-1
					jj = jb+j-1
					kk = kb+k-1
					gstencil(i,j,k,1) = ib+(i-1)
					gstencil(i,j,k,2) = jb+(j-1)   
					gstencil(i,j,k,3) = kb+(k-1)
#if !PARALLEL
					if(ii.lt.1) ii = mxf+ii-1
					if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
					if(jj.lt.1) jj = my+jj
					if(jj.gt.my) jj = jj-my
					if(kk.lt.1) kk = mz+kk
					if(kk.gt.mz) kk = kk-mz

					LOCAL_INDEX(ii)
#if PARALLEL
					if (flag.eq.1) then
						if (ii.lt.0) then 
							PRint*,'fesdata:ERROR WEST SIDE', myid,frombnd,fromrpr
							PRint*,'WEST SIDE: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
							PRint*,'WEST SIDE: INDEX: ', l
							PARALLEL_FINISH()
							STOP
						elseif (ii.gt.nx+1) then
							PRint*,'fesdata:ERROR EAST SIDE', myid,frombnd,fromrpr
							PRint*,'EAST SIDE: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
							PRint*,'EAST SIDE: INDEX: ', l
							PARALLEL_FINISH()
							STOP
						endif
					else
	               if (ii.lt.-1) then 
							PRint*,'fesdata:ERROR WEST EXT POint', myid,frombnd,fromrpr
							PRint*,'WEST EXT POint: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
							PRint*,'WEST EXT POint: INDEX: ', l
							PARALLEL_FINISH()
							STOP
						elseif (ii.gt.nx+2) then
							PRint*,'fesdata:ERROR EAST EXT POint', myid,frombnd,fromrpr
							PRint*,'EAST EXT POint: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
							PRint*,'EAST EXT POint: INDEX: ', l
							PARALLEL_FINISH()
							STOP
						endif
					endif
#endif
					if (uext) then
						vsten(i,j,k,1:ndim) = ur(ii,jj,kk,1:ndim)
					else
						vsten(i,j,k,1:ndim) = urstar(ii,jj,kk,1:ndim)
						ppgrsten(i,j,k,1:ndim) = ppr(ii,jj,kk,1:ndim)
					endif
				enddo
			enddo
		enddo

		call interpolator(gstencil(1:onew,1:onew,1:onew,1:3),vsten(1:onew,1:onew,1:onew,1:ndim),pos(1:ndim),ul(1:ndim),onew,interp_scheme,weightp) 

		if (.not.uext) then
			do n = 1, ndim
				ppll(n) = array_dot_product(ppgrsten(1:onew,1:onew,1:onew,n),weightp(1:onew,1:onew,1:onew))
			end do
		endif
	end subroutine fes_interpolate_data


	subroutine fes_poisson(j,k)
#if PARALLEL
		USE nlarrays, ONLY: westbuf => uf1, eastbuf => uf2
#endif
		implicit none 
		INTEGER, INTENT(IN) :: j,k
		!-----------------------------------------------------------------------
		!	local variables
		real(prcn) ::  a(nx),b(nx),c(nx)

		real(prcn) ::  dv(nx),dvi(nx)
		real(prcn) ::  dv1(nx),dv1i(nx)
		real(prcn) ::  dmax,dmaxi, prftmp, alp_coef_dt
		complex(prcn) ::  tmpc,div(nx), utemp(3), tmpp, tmpc1, tmpc2, tmpc3
		integer ::  i, n, ip, im
		REAL(prcn):: auxu(nx), auxv(nx), alpha, beta, auxz(nx), vdotz,vdotyr,vdotyi, tmpdivr, tmpdivi, rout
		REAL(prcn) :: send_temp(3), recv_temp(3)

		!-----------------------------------------------------------------------
		!	calculate divergence of intermediate velocity field
		!	form divergence in Fourier space on pressure grid

		alp_coef_dt = dt

		alpha = one
		beta = one

		do i=1,nx !mx1
			a(i)=one/dx2-qrt*w2(j,k)
			b(i)=-two/dx2-half*w2(j,k)
			c(i)=one/dx2-qrt*w2(j,k)

			dv(i)=zero
			dvi(i)=zero
			div(i)=czero

			auxu(i) = zero
			auxv(i) = zero
			auxz(i) = zero
		enddo

		do i=1,nx !mx1
			!remember u* is stored in onl
			tmpc=czero
			tmpc=(onl(i+1,j,k,1)-onl(i,j,k,1))/(dx*dt) + (ff(i+1,j,k,1)-ff(i,j,k,1))/dx
			tmpc=tmpc+half*wy(j)*(onl(i,j,k,2)+onl(i+1,j,k,2))/dt + wy(j)*(ff(i+1,j,k,2)+ff(i,j,k,2))/two
			tmpc=tmpc+half*wz(k)*(onl(i,j,k,3)+onl(i+1,j,k,3))/dt + wz(k)*(ff(i+1,j,k,3)+ff(i,j,k,3))/two

			!       multiplying with g(i) to remove source term in poisson equation

			dv(i)= dreal(tmpc)
			dvi(i)= dimag(tmpc)

#if PARALLEL
			! I HAVE NO IDEA WHY. CHECK FOR PARALLEL VERSION
			onl(i,j,k,1) = tmpc/alp_coef_dt
#endif
		enddo

    !-----------------------------------------------------------------------
    !
#if !PARALLEL
		if(xperiodic) then
			if(j.eq.1.and.k.eq.1)then
				call tridag(a(1:nx-1),b(1:nx-1),c(1:nx-1),dv(1:nx-1),dv1(1:nx-1),nx-1)
				dv1(nx) = zero
				dv1i(1:nx) = zero
				vdotz = zero
				vdotyr = zero
				vdotyi = zero
			else
				auxu(1) = alpha
				auxu(nx) = beta
				auxv(1) = c(nx)/auxu(nx)
				auxv(nx) = a(1)/auxu(1)
				b(1) = b(1) - auxu(1)*auxv(1)
				b(nx) = b(nx) - auxu(nx)*auxv(nx)

				call tridag3(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),auxu(1:nx)&
				&,dvi(1:nx),dv1(1:nx),auxz(1:nx),dv1i(1:nx),nx)
				vdotz = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
				vdotyr = auxv(1)*dv1(1) + auxv(nx)*dv1(nx)
				vdotyi = auxv(1)*dv1i(1) + auxv(nx)*dv1i(nx)
			end if
		end if

		do i=1,nx !mx1
			tmpdivr = dv1(i) - vdotyr*auxz(i)/(one + vdotz)
			tmpdivi = dv1i(i) - vdotyi*auxz(i)/(one + vdotz)

			div(i)=dcmplx(tmpdivr,tmpdivi)
		enddo

		!-----------------------------------------------------------------------
		!	correct intermediate velocity field divergence error
		!	
		!	u(t+dt)=u(int)-grad(phi)
		!   u(k) = u* - dt*grad(phi) + dt*f
		!-----------------------------------------------------------------------
    
		if(xperiodic) then
			u(1,j,k,1) = onl(1,j,k,1) - dt*(div(1)-div(mx1))/dx + dt*ff(1,j,k,1)
			u(1,j,k,2) = onl(1,j,k,2) - dt*wy(j)*half*(div(1)+div(mx1)) + dt*ff(1,j,k,2)
			u(1,j,k,3) = onl(1,j,k,3) - dt*wz(k)*half*(div(1)+div(mx1)) + dt*ff(1,j,k,3)
		endif

		do i=2,nx  
			u(i,j,k,1)=onl(i,j,k,1)- dt*(div(i)-div(i-1))/dx + dt*ff(i,j,k,1)
			u(i,j,k,2)=onl(i,j,k,2)- dt*half*wy(j)*(div(i-1)+div(i)) + dt*ff(i,j,k,2)
			u(i,j,k,3)=onl(i,j,k,3)- dt*half*wz(k)*(div(i-1)+div(i)) + dt*ff(i,j,k,3)
		enddo

		if(xperiodic) u(mx,j,k,:) = u(1,j,k,:)

		!-----------------------------------------------------------------------
		!	update pressure field
		!	A positive sign should be here...
		!	P = P* + phi - half*dt*vis*(grad^2 (phi)) + half*dt*vis* grad(f)
		!   Store pressure gradient in ff

		do i=1,nx !mx1!=mx-1
			ip = i+1
			im = i-1

			! CHECK THESE VALUES FOR THE PARALLEL VERSION
			if (i==1) im = mx1
			if (i==nx) ip = 1

!			if(i.eq.1)then
!				p(i,j,k)=pstar(i,j,k)+div(i) - (half*dt)*vis*((div(mx1)-2*div(i)+div(i+1))/dx2 -w2(j,k)*(div(mx1)+2*div(i)+div(i+1))/4.)
!
!				ff(i,j,k,1) = (div(i+1)-div(nx))/(dx) 
!				ff(i,j,k,2) = wy(j)*half*(div(1)+div(nx))
!				ff(i,j,k,3) = wz(k)*half*(div(1)+div(nx))
!			elseif(i.eq.nx)then
!				p(i,j,k)=pstar(i,j,k)+div(i) - (half*dt)*vis*((div(i-1)-2*div(i)+div(1))/dx2 -w2(j,k)*(div(i-1)+2*div(i)+div(1))&
!				&/4.)
!				ff(i,j,k,1) = (div(i)-div(i-1))/(dx) 
!				ff(i,j,k,2) = wy(j)*half*(div(i)+div(i-1))
!				ff(i,j,k,3) = wz(k)*half*(div(i)+div(i-1))
!			else
				! WHY NOT JUST - W2*DIV(I)
!				p(i,j,k)=pstar(i,j,k)+div(i) - half*dt*vis * ( (div(im)-2*div(i)+div(ip))/dx2 - w2(j,k) * (div(im)+2*div(i)+div(ip))/4.)
				p(i,j,k)=pstar(i,j,k)+div(i) - (half*dt)*vis*((div(im)-2*div(i)+div(ip))/dx2 - w2(j,k)*div(i))
				p(i,j,k) = p(i,j,k) + half*dt*rhof*vis * ((ff(i+1,j,k,1)-ff(i,j,k,1))/dx + half*wy(j)*(ff(i+1,j,k,2)+ff(i,j,k,2)) + half*wz(k)*(ff(i+1,j,k,3)+ff(i,j,k,3)))

				ff(i,j,k,1) = (div(i)-div(i-1))/(dx) 
				ff(i,j,k,2) = wy(j)*half*(div(i)+div(i-1))
				ff(i,j,k,3) = wz(k)*half*(div(i)+div(i-1))
!			endif
		enddo
		! CHECK THESE VALUES FOR THE PARALLEL VERSION
		ff(mx,j,k,:) = ff(1,j,k,:)
#endif
		RETURN
	end subroutine fes_poisson


	SUBROUTINE correct_mean_quantities
		USE nlmainarrays, ONLY : ur=>ubcp
		IMPLICIT NONE
		REAL(prcn) :: mean_force(ndim), phase_mass, slip_obt(ndim) 
		INTEGER :: idim, m

		umean(:) = umean_star(:) + dt*frmean(:)
		phase_mass = rhos*pi*dia_phys**3.d0/6.d0
    

		!SHOULD FORCE BE PRES+VISC OR PRES+VISC-MPG*V
		do idim = 1, ndim
			mean_force(idim) = zero
			do m = 1,nbody 
				mean_force(idim) = mean_force(idim) + force(m,idim)
			enddo 
			mean_force(idim) = mean_force(idim)/real(nbody,prcn)

!			if (move_particles) usmean(idim) = VMEAN_STAR(idim) + half*DT/phase_mass*mean_force(idim)
		enddo

		CALL calc_velreal(u, umean, ur)
		CALL compute_mean_fluid_velocity
#if 0
		WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN (k) = ', (UFMEAN(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'USMEAN (k) = ', (USMEAN(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'MEAN FORCE (k) = ', (MEAN_FORCE(idim), idim = 1, ndim)
#endif

		if(TRIM(simtype).eq.'sediment')then
			afk(:) = zero
			betak(:) = rhof/DT*(ufmean(:)-ufmean_des(:))
		else if(TRIM(simtype).eq.'riser')then
			afk(:) = zero
			betak(:) = zero   
		else if(TRIM(simtype).eq.'accframe')then
			slip_obt(:) = usmean(:) - ufmean(:)
			slip_des(:) = usmean_des(:) - ufmean_des(:)

			betak(:) = (slip_obt(:) - slip_des(:))/DT
			betak(:) = betak(:)/(one/rhos - one/rhof) 

			afk(:) = (usmean(:)-usmean_des(:))/DT- betak(:)/rhos
		else if(TRIM(simtype).eq.'fixed')then
			afk(:) = zero

			slip_obt(:) = usmean(:) - ufmean(:)
			slip_des(:) = usmean_des(:) - ufmean_des(:)

			betak(:) = - (slip_obt(:) - slip_des(:)) * rhof/dt
!			betak(:) = betak(:)/(one/rhos - one/rhof) 
		end if
    
		umean(:) = umean_star(:) + dt * (-betak(:)/rhof - afk(:) + frmean(:))

#if 0
		WRITE(*,'(A25,3(2x,g17.8))')'BETA (k) = ', (BETAK(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'AF (k) = ', (AFK(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'<f> (k) = ', (FRMEAN(idim), idim = 1, ndim)
		WRITE(*,'(A25,3(2x,g17.8))')'<U>(k) = ', (UMEAN(idim), idim = 1, ndim)
#endif

	END SUBROUTINE correct_mean_quantities


	SUBROUTINE CHECK_ITER_CONVERGENCE(iterstop)
		IMPLICIT NONE
		LOGICAL, Intent(out) :: iterstop
		INTEGER :: idim
		REAL(prcn) :: mean_force(ndim), forcemod, ferror
		real(prcn), save :: tmp0, tmp1, tmp2, beta_error

		if (itercount==1) then
			tmp0 = sqrt(dot_product(betak,betak))
			tmp1 = zero
		else
			tmp1 = tmp2
		endif
		tmp2 = sqrt(dot_product(betak,betak))
		beta_error = abs(tmp2-tmp1)/tmp0

		if (move_particles) then
			do idim = 1, ndim
				mean_force(idim) = SUM(force(1:nbody,idim))/real(nbody,prcn)
			end do

			forcemod = DOT_PRODUCT(mean_force(1:ndim), mean_force(1:ndim))
			forcemod = DSQRT(forcemod)

			if(forcemod.gt.zero)then
				ferror = ABS(forcemod - forcemod_old)/forcemod
			else
				ferror = one
			end if
			forcemod_old = forcemod
			Write(*,'(A,2x,I,3(2x,d15.7))')'RESIDUAL IN INNER ITERATIONS: ', itercount, ferror, beta_error ,forcemod/(3.d0*pi*vis*dia_phys*ucharmod)

			iterstop = (ferror.lt.1E-05).or.(itercount.gt.INNERITERMAX)
		else
			write(*,'(A,2x,I,3(2x,d15.7))')'RESIDUAL IN INNER ITERATIONS: ', itercount, beta_error

			iterstop = (beta_error.lt.1E-03).or.(itercount.gt.INNERITERMAX)
		endif
		!iterstop = (itercount.eq.3)
		iterstop = iterstop.and.(itercount.ge.1)
	END SUBROUTINE CHECK_ITER_CONVERGENCE


	SUBROUTINE FES_COMPUTE_NEW_TIMESTEP(check_dt)
		USE nlmainarrays, Only : ubcp, pbcp
		USE nlarrays, Only : uf1, uf2, uf3
		USE scalar_data
		IMPLICIT NONE

		logical, intent(in) :: check_dt

		REAL(PRCN) ::  u_max, v_max, w_max, umax_tmp
		REAL(PRCN) ::  umax_loc, vmax_loc, wmax_loc, mixmeanslip(ndim), umean_temp(ndim)
		INTEGER :: idim, i, j, k, iphs, partstart, partend

		if (.not.check_dt) goto 10

		DTNEW = LARGE_NUMBER
		u_max = SMALL_NUMBER
		v_max = SMALL_NUMBER
		w_max = SMALL_NUMBER
		umax_loc = SMALL_NUMBER
		vmax_loc = SMALL_NUMBER
		wmax_loc = SMALL_NUMBER

		umax_tmp = LARGE_NUMBER
       
#if PARALLEL
		do i = 1, nx
#else
		do i = 1,nx
#endif
			do k = 1, mz
				do j = 1, my 
					if((i.gt.0).and.(i.lt.nx+1))then
						umax_loc = MAX(umax_loc, ABS(ubcp(i,j,k,1)-mesh_vel(1)))
						vmax_loc = MAX(vmax_loc, ABS(ubcp(i,j,k,2)-mesh_vel(2)))
						wmax_loc = MAX(wmax_loc, ABS(ubcp(i,j,k,3)-mesh_vel(3)))
					endif
				enddo
			enddo
		enddo
		GLOBAL_DOUBLE_MAX(umax_loc,u_max,1,decomp_group)
		GLOBAL_DOUBLE_MAX(vmax_loc,v_max,1,decomp_group)
		GLOBAL_DOUBLE_MAX(wmax_loc,w_max,1,decomp_group)
       
		umax_tmp = MAX(u_max,v_max)
		umax_tmp = MAX(umax_tmp,w_max)

		WRITE(*,'(A25,2(2x,g12.5))') "MAX VEL AND CFL:", umax_tmp,  umax_tmp*dt/dx
		!Write(*,*)'umax: ', umax_tmp
		if(umax_tmp*dt/dx.gt.cfl) then
			dtnew = cfl*dx/umax_tmp
!			dtchanged = .TRUE.
		else
			dtnew = dt 
!			dtchanged = .FALSE.
		endif

		IF(adaptive.and..not.only_dem)then
			DT = DTNEW
		ENDIF
       cf = -one/dt
       cforig = cf
       
		if(I_AM_NODE_ZERO) WRITE(*,'(A, i4,2(2x,g12.5),/,2(A30, g12.5,/))')'IDUMSTEP, t, tend = ', IDUMSTEP, t, tendused ,& 
			'cf original = ', cforig, &
			'cf used = ', cf


		!usmean can become anything in IBM. SO we dont evaluate mean slip based on the current value of usmean. As long as the particle velocities are correctly attained, usmean_des has been reached and usmean calculated is irrelevant.!
10		continue
		usmean = zero
		do iphs = 1, nphases
			phase_array(iphs)%mean_spec_vel(1:ndim) = zero
			partstart = phase_array(iphs)%pstart
			partend = phase_array(iphs)%pend
			do idim = 1, ndim
				phase_array(iphs)%mean_spec_vel(idim) = SUM(velbdy(partstart:partend,idim))/real(phase_array(iphs)%npart,prcn)
			end do
!!$		usmean(1:ndim) = usmean(1:ndim) + (phase_array(iphs)%volfracg)*(phase_array(iphs)%mean_spec_vel(1:ndim))
			usmean(1:ndim) = usmean(1:ndim) + (phase_array(iphs)%volfrac)*(phase_array(iphs)%mean_spec_vel(1:ndim))
			!WRite(*,*)'volfracg = ', (phase_array(iphs)%volfracg)
		end do

		usmean(:) = usmean(:)/mean_volfrac

		meanslip(:) = (one-maxvolfrac)*(usmean_des(:)-ufmean(:))
		meanslipmod = DSQRT(dot_product(meanslip(1:ndim), meanslip(1:ndim)))

		mixmeanslip(:) = (one-maxvolfrac)*(usmean(:)-ufmean(:))
		mixmeanslipmod = DSQRT(DOT_PRODUCT(mixmeanslip(1:ndim),mixmeanslip(1:ndim)))

		if (I_AM_NODE_ZERO) then
			WRITE(*,'(A25,3(2x,g17.8))')'USMEAN DES = ', usmean_des(:)
			WRITE(*,'(A25,3(2x,g17.8))')'USMEAN ACT = ', usmean_act(:)
			WRITE(*,'(A25,3(2x,g17.8))')'USMEAN MIX = ', usmean(:)
			WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN DES = ', ufmean_des(:)
			WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN ACTUAL = ', ufmean(:)
			!WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN CC = ', ufmean_cc(:)
			WRITE(*,'(A25,3(2x,g17.8))')'UMEAN ACTUAL = ', umean(:)
			WRITE(*,'(A25,2(2x,g12.5))') "MAX VEL AND CFL:", umax_tmp,  umax_tmp*dt/dx
			WRITE(*,'(A40,3(2x,g12.5))') "MEANSLIPMOD AND REY(MEANSLIPMOD):", mixmeanslipmod, mixmeanslipmod*char_length/vis
		endif
	END SUBROUTINE FES_COMPUTE_NEW_TIMESTEP
end module fes
