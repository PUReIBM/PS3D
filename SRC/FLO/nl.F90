Module nlcalc 
#include "ibm.h"
	use precision 
	use constants 
	use global_data
	use machine
	Use boundary_condition
	Use dependent_functions
	Use collision_mod
	use fftw3_interface
	implicit none 
contains   
	subroutine nonlinear(rks)
		use bcsetarrays, only : fr
		use nlmainarrays, Only :  ubc, nlbc, onlbc
		use init_turb, only : forced_turbulence, linear_forcing, eddy_time_i, meaniso_force, activate_forcing_time

		implicit none 
		integer, Intent(in) :: rks
		real(prcn) :: nlmean(ndim), nlmean_loc(ndim) 

		integer i,j,k, idim
		real(prcn) :: utemp, tempt, dt_fluid_orig

		real(prcn) cpu0, cpu1, frmeanloc(ndim)

		if (I_AM_NODE_ZERO) call CPU_TIME (CPU0)

		!-----------------------------------------------------------------------
		!	reset terms to be overwritten
		onlbc(:,:,:,:) = nlbc(:,:,:,:)
		onl(:,:,:,:) = nl(:,:,:,:)
		nl(:,:,:,:) = czero

		if (I_AM_NODE_ZERO) write(*,'(A50)')'CONSERVATIVE: performing del dot uu'
		if (.not.only_dem) call form_nl

		t=t+dt*(coef(rks,1)+coef(rks,2))

		if (move_particles) S_TIME = t

		utemp = ucharmod 

		tempt = ramp_frac_time*float(global_n(1))*dx/utemp

		if (move_particles) soft_start=.false.

		if (soft_start) then 
			if (t.LT.tempt) then
				cf = -one/dt
				cforig = cf
				cf=cf*(t/tempt)
			else
				cf = -one/dt
				cforig=cf
			endif
		else
			cf = -one/dt
			cforig = cf
		endif
    
		if (I_AM_NODE_ZERO) write (*,'(A, i4,2(2x,g12.5),/,2(A30, g12.5,/))') 'IDUMSTEP, t, tend = ', &
									& IDUMSTEP, t, tendused , 'cf original = ', cforig, 'cf used = ', cf
    
		!Jamals convention of absorbing the negativesign into the nl term 
		if (aliasflag.eq.1) then
			!new nl(1:nx+1,:,:,n) = -thrd*nl(1:nx+1,:,:,n) !Exposing the last stride
			nl(:,:,:,:) = -thrd*nl(:,:,:,:) 
		else
			!new nl(1:nx+1,:,:,n) = -nl(1:nx+1,:,:,n)
			nl(:,:,:,:) = -nl(:,:,:,:)
		endif
    
		!-----------------------------------------------------------------------
		!	store nonlinear terms in forcing domain
		do idim=1, ndim
			call fftwc2r( nl(:,:,:,idim),  nlbc(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))
			call communicate_in_gohst_domain( nlbc(0:local_ni(1)+1,0:local_ni(2)+1,0:local_ni(3)+1,idim))

			if (iglobstep==1) then
				call fftwc2r(onl(:,:,:,idim), onlbc(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))
				call communicate_in_gohst_domain(onlbc(:,:,:,idim))
			endif
		enddo

		if (debug_check) then
			nlmean(1:ndim) = zero
			do idim = 1, ndim
				nlmean_loc(idim) = SUM(nlbc(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))
				GLOBAL_DOUBLE_SUM(nlmean_loc(idim),nlmean(idim),1,comm_cart_2d)
			enddo

			nlmean(:) = nlmean(:)/(global_n(1)*global_n(2)*global_n(3))
			if (I_AM_NODE_ZERO) write (*,'(A25,3(2x,g12.5))') 'NLMEAN = ', nlmean
		endif
    
		if (I_AM_NODE_ZERO) then
			call CPU_TIME (CPU1) 
			nl_time_current = cpu1-cpu0
			nl_time = nl_time + nl_time_current
		endif

		if (nbody.gt.0) then
			call bcset(rks)
		endif

		if (trim(input_type).eq."single-phase") then
			if (iturbon.and.forced_turbulence) then
				if (t>activate_forcing_time*eddy_time_i) then
					if (I_AM_NODE_ZERO) write (*,*) "LINEAR_FORCING ACTIVE...."
					call linear_forcing

					do idim = 1, ndim 
						fr(:,:,:,idim) = fr(:,:,:,idim) - meaniso_force(idim)
						call fftwr2c(fr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim), ff(:,:,:,idim))
					enddo
				endif
			endif
		endif

		if (I_AM_NODE_ZERO) write (*,'(A)') 'END OF NONLINEAR'
	end subroutine nonlinear

	subroutine form_nl
		use nlmainarrays, Only : ubcp
		use field_tmp_arrays
		implicit none 
		integer :: i,j,k, dim1, dim2, ii, jj, kk
		real(prcn) :: slope_factor, mesh_veltemp(ndim)
		real(prcn) ::  u_max, v_max, w_max
		real(prcn) ::  umax_loc, vmax_loc, wmax_loc
		complex(prcn) :: wtmp1, wtmp2

		slope_factor = one !nl_slope_factor
		!phi_ri = one
		if (I_AM_NODE_ZERO) write (*,'(A,2x,g17.8)') 'SLOPE FACTOR IN NL = ', slope_factor 

		mesh_veltemp = mesh_vel
		if (.not.movingcv) mesh_vel = zero
		if (I_AM_NODE_ZERO) write (*,'(A,3(2x,g17.8))')'MESH VEL IN NL = ', mesh_vel(:) 

		do dim1=1, ndim
			do dim2=1, dim1
				do k=1, local_ni(3)
					do j=1, local_ni(2)
						do i = 1, local_ni(1)
							urtmp(i,j,k) = (ubcp(i,j,k,dim1)-mesh_vel(dim1))*(ubcp(i,j,k,dim2)-mesh_vel(dim2))
						enddo
					enddo
				enddo

				call fftwr2c(urtmp(1:local_ni(1),1:local_ni(2),1:local_ni(3)), uftmp)

#if !PARALLEL
				do k=1, local_no(3)
					do j=1, local_no(2)
						do i=1, local_no(1)
							if (dim1==1) then
								wtmp1 = wx(i)
							elseif (dim1==2) then
								wtmp1 = wy(j)
							elseif (dim1==3) then
								wtmp1 = wz(k)
							endif

							if (dim2==1) then
								wtmp2 = wx(i)
							elseif (dim2==2) then
								wtmp2 = wy(j)
							elseif (dim2==3) then
								wtmp2 = wz(k)
							endif
#else
				!REMEMBER DATA IS TRANSPOSED IN FOURIER SPACE
				do k=1, local_no(2)
					do j=1, local_no(1)
						do i=1, local_no(3)
							if (dim1==1) then
								wtmp1 = wy(j) !DERIVATIVE ALONG X
							elseif (dim1==2) then
								wtmp1 = wz(k) !DERIVATIVE ALONG Y
							elseif (dim1==3) then
								wtmp1 = wx(i) !DERIVATIVE ALONG Z
							endif

							if (dim2==1) then
								wtmp2 = wy(j) !DERIVATIVE ALONG X
							elseif (dim2==2) then
								wtmp2 = wz(k) !DERIVATIVE ALONG Y
							elseif (dim2==3) then
								wtmp2 = wx(i) !DERIVATIVE ALONG Z
							endif
#endif
							nl(i,j,k,dim1) = nl(i,j,k,dim1) + uftmp(i,j,k) * wtmp2
							if (dim1/=dim2) nl(i,j,k,dim2) = nl(i,j,k,dim2) + uftmp(i,j,k) * wtmp1
						enddo
					enddo
				enddo
			enddo
		enddo
		mesh_vel = mesh_veltemp
	end subroutine form_nl
end Module nlcalc













