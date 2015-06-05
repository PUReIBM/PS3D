MOdule mypost_process
#include "../FLO/ibm.h"
	use precision
	use global_data
	use nlmainarrays , only : ubcp, pbcp
	use fftw3_interface
	use constants
	use randomno
	use general_funcs
	use string_funcs
	!use boundary_condition
	!use dependent_functions
	use bcsetarrays
	use field_tmp_arrays
!	use epsarray
	use postproc_funcs
!	use string_funcs
!	use init_turb
	use parallel
	implicit none
!	private
!	public :: dissip_continuous, compute_sijsij

	real(prcn) :: ufmean_mag, usmean1_mag, usmean2_mag, mix_mean_vel_mag, mix_mean_slip_mag, sijsij, sijsij_pure, tke_pure, source, dissipation, int_tke_pure !tke

	real(prcn), allocatable :: mix_mean_slip(:), mean_spec_vel_mag(:), &
		& mix_mean_vel(:), grant(:), tke_s_pure(:), acc_mean(:), acc_variance(:), vel_mean(:), uf_cage(:,:)

!	integer :: nvar, nvar1, nvar2, 
	integer :: line=1

	logical :: first_pass_post

!	type :: slice_type
!		integer :: num, dy
!		integer, allocatable :: ind(:), fluidnum(:), solidnum(:)
!		real(prcn), allocatable :: ufmean(:,:), mpg(:,:)
!		real(prcn), allocatable :: areaf(:)
!	end type slice_type
!	type(slice_type),allocatable:: slice(:)

	real(prcn), allocatable :: p_acc_avg(:)
	integer :: p_acc_avg_count

contains
#if 0
	subroutine compute_mix_mean_vel
		implicit none
		integer :: pstart, pend, m, iphs, idim

		if (.not.allocated(mix_mean_vel)) then
			allocate(mix_mean_vel(ndim))
			allocate(mean_spec_vel_mag(nphases), mix_mean_slip(ndim))
			allocate(grant(nphases))
		endif

		mix_mean_vel  = zero
		grant = zero

		do iphs = 1, nphases
			pstart = phase_array(iphs)%pstart
			pend = phase_array(iphs)%pend
			mean_spec_vel_mag(iphs) = sqrt(dot_product(phase_array(iphs)%mean_spec_vel(:),&
										& phase_array(iphs)%mean_spec_vel(:)))

			mix_mean_vel(:) =  mix_mean_vel(:) + phase_array(iphs)%mean_spec_vel(:) * phase_array(iphs)%volfracg / mean_volfrac


			do m = pstart, pend
				do idim = 1, ndim
					grant(iphs) = grant(iphs) + (phase_array(iphs)%mean_spec_vel(idim)-velbdy(m,idim))**2.d0
				end do
			end do
			grant(iphs) = grant(iphs)/(three*phase_array(iphs)%npart)
		end do

		mix_mean_vel(:) = usmean(:)

		mix_mean_slip(:) = mix_mean_vel(:)-ufmean(:)
		mix_mean_vel_mag = dot_product(mix_mean_vel,mix_mean_vel)
		mix_mean_slip_mag = sqrt(dot_product(mix_mean_slip(:), mix_mean_slip(:)))
		pstart = pend+1
	end subroutine compute_mix_mean_vel
#endif


	SUBROUTINE calc_anisotropy(uij, bij, xi, eta)
		implicit none

		REAL(prcn),Intent(in)  :: uij(ndim,ndim)
		REAL(prcn),Intent(out) :: bij(ndim,ndim), xi, eta
		Integer :: m,i,j,k
		Real(prcn) :: trace, aij(ndim,ndim), delta_ij, tmp

		trace = uij(1,1)+uij(2,2)+uij(3,3)
		!PRINT*,' trace =', trace 
		if(.not.(ABS(trace).GT.SMALL_NUMBER)) then
			PRINT*,'TRACE OF THE TENSOR IS VERY SMALL', trace
			RETURN
		end if
		do i = 1, ndim
			do j = 1, ndim
				if(i.eq.j)then 
					delta_ij = 1
				else
					delta_ij = 0
				endif
				bij(i,j) = uij(i,j)/trace - delta_ij/3.d0
			enddo
		enddo

		eta = zero
		xi = zero

		do i = 1,ndim
			do j = 1, ndim
				eta = eta + bij(I,J)*bij(j,i)
				do k= 1, ndim
					xi = xi + bij(i,j)*bij(j,k)*bij(k,i)
				end do
			end do
		end do

		eta = sqrt(eta/6.d0)
		xi = xi/6.d0

		tmp = abs(xi)**(1./3.)
		xi = sign(tmp, xi)
	end subroutine calc_anisotropy

	subroutine reynolds_stress_tensor
		implicit none
		real(prcn) :: spec_fluc(nphases), uiuj_s(nphases,ndim,ndim), tke_s(nphases), tke_mix_s
		real(prcn), dimension(ndim,ndim) :: uiuj_f, bij_f, bij_s
		real(prcn) :: xi, eta1, norm

		real(prcn) :: tmp_tensor(ndim,ndim) !, tmp_vec(ndim)
		integer :: dim1, dim2
		integer :: i, j, k, idim, pstart, pend, m, iphs
		character*100 filename
		integer :: unitnum

		integer :: nvar, nvar1, nvar2
		logical :: filexist

		if (I_AM_NODE_ZERO) write (*,*) "IN REYNOLDS_STRESS_TENSOR"
		if (.not.post_no_flow_mem_alloc) then
			!^^^^^^^ Computing RSM_F ^^^^^^^^^^^^^^
			uiuj_f = 0d0
			do k=1, local_ni(3)
				do j=1, local_ni(2)
					do i=1, local_ni(1)
						if (fluid_atijk(i,j,k)) then
							do dim1=1, ndim
								do dim2=1, dim1
									uiuj_f(dim1,dim2) = uiuj_f(dim1,dim2) + (ubcp(i,j,k,dim1)-ufmean(dim1)) &
																					 *  (ubcp(i,j,k,dim2)-ufmean(dim2))
									if (dim2/=dim1) then
										uiuj_f(dim2,dim1) = uiuj_f(dim2,dim1) + (ubcp(i,j,k,dim1)-ufmean(dim1)) &
																						 *  (ubcp(i,j,k,dim2)-ufmean(dim2))
									endif
								enddo
							enddo
						endif
					enddo
				enddo
			enddo

			!^^^^^^^ Computing RSM_S ^^^^^^^^^^^^^^
			spec_fluc = zero
			uiuj_s    = zero
			tke_s     = zero
			tke_mix_s = zero

			if (move_particles) then
				pstart = 1
				do iphs = 1, nphases
					pend = pstart + phase_array(iphs)%npart - 1
					do m = pstart, pend
						do dim1 = 1, ndim
							do dim2=1, dim1
								uiuj_s(iphs,dim1,dim2) = uiuj_s(iphs,dim1,dim2) + &
											&	(velbdy(m,dim1)-phase_array(iphs)%mean_spec_vel(dim1)) * &
											&	(velbdy(m,dim2)-phase_array(iphs)%mean_spec_vel(dim2))
								if (dim1/=dim2) then
									uiuj_s(iphs,dim2,dim1) = uiuj_s(iphs,dim2,dim1) + &
												&	(velbdy(m,dim1)-phase_array(iphs)%mean_spec_vel(dim1)) * &
												&	(velbdy(m,dim2)-phase_array(iphs)%mean_spec_vel(dim2))
								endif
							enddo
						enddo
					enddo
					uiuj_s(iphs,:,:) = uiuj_s(iphs,:,:)/real(phase_array(iphs)%npart,prcn)
					pstart = pend + 1
				enddo

				do iphs = 1, nphases
					do idim=1, ndim
						tke_s(iphs) = tke_s(iphs) + uiuj_s(iphs,idim,idim)
					enddo
				enddo

				tke_s(:) = tke_s(:)/2
			endif

			if (.not.allocated(tke_s_pure)) allocate(tke_s_pure(nphases))

			tke_s_pure = tke_s

			do iphs = 1, nphases
				tke_mix_s = tke_mix_s + tke_s(iphs)*phase_array(iphs)%volfracg/mean_volfrac
			end do

			tmp_tensor = uiuj_f
			GLOBAL_DOUBLE_SUM(tmp_tensor,uiuj_f,9,comm_cart_2d)

			uiuj_f   = uiuj_f / count_fluid
			tke      = (uiuj_f(1,1)+uiuj_f(2,2)+uiuj_f(3,3)) / 2
			tke_pure = tke

			if (zero_slip) then
				if (iturbon) then
					norm = tke_i
				elseif (ReT>small_number) then
					norm = 3*grant_i/2
				endif
			else
				norm = half * dot_product(usmean_des(:)-ufmean(:),usmean(:)-ufmean(:))
			endif


			if(I_AM_NODE_ZERO)then
				call screen_separator(40,'K')
				do iphs = 1, nphases
					call screen_separator(20,'^')
					write (*,"(1A,1I3)") "SPECIES #", iphs
					!write (*,"(1A,4D15.7)") "MEAN VELOCITY, MAG = ", phase_array(iphs)%mean_spec_vel(:), sqrt(dot_product(phase_array(iphs)%mean_spec_vel(:), phase_array(iphs)%mean_spec_vel(:)))
					write (*,"(1A)") "REYNOLDS STRESS"
					do idim=1, ndim
						write (*,"(3D15.7)") uiuj_s(iphs,idim,:)/norm
					enddo
					!call calc_anisotropy(uiuj_s, bij_s, xi, eta1)
					write (*,"(1A,1D15.7)") "TKE_S = ", tke_s(iphs)/norm
					call screen_separator(20,'-')
				enddo

				write (*,"(1A,1D15.7)") "MIX TKE SOLID = ", tke_mix_s/norm

				call calc_anisotropy(uiuj_f, bij_f, xi, eta1) 
				write (*,'(A,3D15.7)') "TKE_F/MEAN_energy,XI, ETA = ", tke/norm, xi, eta1
				write (*,'(A)') "ANISOTROPY TENSOR = "
				write (*,'(3D15.7)') ((bij_f(i,j), j=1,ndim), i=1,ndim)
				call screen_separator(40,'K')

				filename = trim(run_name)//"_tke"
				unitnum = 1
				call instant_file_opener(filename,unitnum,.true.)

				if (from_post) then
					if (.not.imove==1) then
						if (nphases==1) then
							write (unitnum,"(4F10.4,21D15.7)") dbydx, lybyd, mean_volfrac, re_out,	&
							& t/t_conv, tke, tke/norm, &
							& ((uiuj_f(i,j)/norm,i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim)
						else
							write (unitnum,"(6f10.4,21D15.7)") mean_volfrac, &
							& (phase_array(iphs)%volfracg, iphs=1,nphases), &
							& (yalpha(iphs), iphs=1,nphases), re_out,	&
							& t/t_conv, tke, tke/norm, &
							& ((uiuj_f(i,j)/norm,i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim)
						endif
					else
						if (nphases==1) then
							write (unitnum,"(6f10.4,31D15.7)") dbydx, lybyd, mean_volfrac, &
							& re_out, rhos, coeff_rest, &
							& t/t_conv, tke,tke/norm, &
							& ((uiuj_f(i,j)/norm,i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim), &
							& tke_s(1)/norm, &
							& ((uiuj_s(1,i,j)/norm,i=1,ndim),j=1,ndim)
						else
							write (unitnum,"(8f10.4,42D15.7)") mean_volfrac,(phase_array(iphs)%volfracg, iphs=1,nphases), (yalpha(iphs), iphs=1,nphases), &
							& re_out, rhos, coeff_rest, &
							& t/t_conv, tke, tke/norm, &
							& ((uiuj_f(i,j)/norm,i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim), &
							& tke_s(1)/norm, &
							& ((uiuj_s(1,i,j)/norm,i=1,ndim),j=1,ndim), &
							& tke_s(2)/norm, &
							& ((uiuj_s(2,i,j)/norm,i=1,ndim),j=1,ndim), &
							& tke_mix_s/norm
						endif
					endif
				else
					write (unitnum,"(17d15.7)") t/t_conv, tke, tke/norm, xi, eta1, &
					&		uiuj_f(1,1), uiuj_f(2,2), uiuj_f(3,3), &
					&		uiuj_f(1,2), uiuj_f(1,3), uiuj_f(2,3), &
					&		bij_f(1,1), bij_f(2,2), bij_f(3,3), &
					&		bij_f(1,2), bij_f(1,3), bij_f(2,3)
				endif
				close (unitnum)
#if 0
				if (imove==1) then
					filename = trim(run_name)//"_tke_mix.dat"
					open (unit=1,file=trim(filename),status="replace")
					write (1,"(4f8.4,4D15.7)") mean_volfrac, re_out, rhos, coeff_rest, tke, tke_mix_s, (tke*(1-mean_volfrac)*rhof + tke_mix_s*mean_volfrac*rhos), (tke*(1-mean_volfrac)*rhof + tke_s*mean_volfrac*rhos) / ((1-mean_volfrac)*rhof + mean_volfrac*rhos)
					close (1)
				endif
#endif
			endif
		else
			if(I_AM_NODE_ZERO)then
				if (.not.imove==1) then
					if (nphases==1) then
						nvar = 25
						nvar1 = 4
						nvar2 = 21
					else
						nvar = 27
						nvar1 = 6
						nvar2 = 21
					endif
				else
					if (nphases==1) then
						nvar = 37
						nvar1 = 6
						nvar2 = 31
					else
						nvar = 50
						nvar1 = 8
						nvar2 = 42
					endif
				endif
				
				filename="_tke_post.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line, include_nmis)
			endif 
		endif
	end subroutine reynolds_stress_tensor



	subroutine compute_sijsij
		use init_turb, only : epsf_forcing, forced_turbulence, eddy_time_i, sampling_dissip_time
		implicit none
	
		real(prcn) :: sijsij_loc
		character*100 filename
		integer :: unitnum
		
	
		real(prcn) :: norm
		integer :: i, j, k, idim
		integer :: dim1, dim2

		integer :: nvar, nvar1, nvar2
		complex(prcn) :: wtmp1, wtmp2
		logical :: filexist

		if (.not.post_no_flow_mem_alloc) then
			!call compute_mix_mean_vel

			sijsij_loc = zero
			if (I_AM_NODE_ZERO) write (*,*) "GENERATING SijSij"

			do dim1=1, ndim
				do dim2=1, dim1
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

								uftmp(i,j,k) = half * (u(i,j,k,dim1)*wtmp2 + u(i,j,k,dim2)*wtmp1)
							enddo
						enddo						
					enddo

					call fftwc2r(uftmp,urtmp)

					do k=1, local_ni(3)
						do j=1, local_ni(2)
							do i=1, local_ni(1)
								if (fluid_atijk(i,j,k)) then
									if (dim1==dim2) then
										sijsij_loc = sijsij_loc + 2*vis *   urtmp(i,j,k)**2
									else
										sijsij_loc = sijsij_loc + 2*vis * 2*urtmp(i,j,k)**2
									endif
								endif
							enddo
						enddo
					enddo
				enddo
			enddo

			sijsij = zero
			GLOBAL_DOUBLE_SUM(sijsij_loc,sijsij,1,comm_cart_2d)

			sijsij_pure = sijsij / global_n(1)/global_n(2)/global_n(3)
			sijsij      = sijsij / global_n(1)/global_n(2)/global_n(3)

			if (iturbon.and.forced_turbulence) then
				if (t>sampling_dissip_time*eddy_time_i .and. epsf_forcing<small_number) then
					!epsf_forcing = sijsij_pure
					epsf_forcing = epsf_i
					if (I_AM_NODE_ZERO) then
						filename = trim(run_name)//"_linear_forcing.dat"
						open (unit=1, file=trim(filename), status="replace", action="write")
						write (1,*) epsf_forcing
						close(1)
					endif
				endif
			endif

			if (I_AM_NODE_ZERO) then
				filename=trim(run_name)//"_sijsij"
				unitnum = 1
				call instant_file_opener(filename,unitnum,.true.)

				if (zero_slip) then
					if (iturbon) then
						norm = epsf_i
					elseif (ReT>small_number) then
						norm = vis * (ucharmod / dia_phys)**2
					endif
				else
					norm = vis * dot_product(usmean(:)-ufmean(:),usmean(:)-ufmean(:)) / dia_phys**2
				endif

				if (from_post) then
					write (unitnum,"(4f14.8,9D15.7)") dbydx, lybyd, re_out, mean_volfrac, &
					& t/t_conv, sijsij, sijsij/norm
				else
					write (unitnum,"(1f14.8,9D15.7)") t/t_conv, sijsij, sijsij/norm
				endif
				close (unitnum)
				write (*,"(1a,1f14.8,9D15.7)") "DISSIPATION = ", t/t_conv, sijsij, sijsij/norm
				call screen_separator(40,'E')
			endif
		else
			if(I_AM_NODE_ZERO)then
				if (.not.imove==1) then
					if (nphases==1) then
						nvar  = 7
						nvar1 = 4
						nvar2 = 3
					else
						nvar  = 21
						nvar1 = 8
						nvar2 = 13
					endif
				else
					if (nphases==1) then
						nvar  = 15
						nvar1 = 6
						nvar2 = 9
					else
						nvar  = 23
						nvar1 = 10
						nvar2 = 13
					endif
				endif
				filename = "_sijsij_post.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line, include_nmis)

				if (.not.imove==1) then
				!	if (nphases==1) then
				!		nvar  = 6
				!		nvar1 = 4
				!		nvar2 = 2
				!	else
				!		nvar  = 14
				!		nvar1 = 8
				!		nvar2 = 6
				!	endif
				!	filename = "_force_comp_post.dat"
				!	call mis_average(nvar, nvar1, nvar2, filename, line, include_nmis)
				endif
			endif 
		endif	
	end subroutine compute_sijsij

	subroutine write_drag
		implicit none
		character*100 filename1, filename2
		integer :: line, nvar, nvar1, nvar2, iphs, unitnum, unitnum_chem
		real(prcn) :: tmp1, tmp2, tmp3, tmp4, ibm_drag


		if (I_AM_NODE_ZERO) then
			write (*,*) "IN WRITE_DRAG..."
			if (.not. post_no_flow_mem_alloc) then
				filename1=trim(run_name)//"_drag"
				filename2=trim(run_name)//"_drag_chem"

				unitnum = 1
				unitnum_chem = 2
				call instant_file_opener(filename1,unitnum,.true.)
				call instant_file_opener(filename2,unitnum_chem,.true.)

				if (nphases==1) then
					call compute_ibm_drag(mean_volfrac, re_out, ibm_drag)

					write (unitnum,"(5f10.4,4D15.7)") dbydx, lybyd, re_out, archno, mean_volfrac, &
					& 	norm_drag, norm_drag_spec(1), ibm_drag, abs(norm_drag_spec(1)-ibm_drag) / ibm_drag

					write (unitnum_chem,"(5f10.4,2D15.7)") dbydx, lybyd, re_out, archno, mean_volfrac, &
					& 	norm_drag_chem, norm_drag_chem_spec(1)
				else
					write (unitnum,"(9f10.4,3D15.7)") dbydx, lybyd, re_out, archno, mean_volfrac, &
					&	(phase_array(iphs)%volfracg, iphs=1,nphases), &
					&	(yalpha(iphs),iphs=1,nphases), 			&
					&	norm_drag, norm_drag_spec(1:nphases)

					write (unitnum_chem,"(9f10.4,3D15.7)") dbydx, lybyd, re_out, archno, mean_volfrac, &
					&	(phase_array(iphs)%volfracg, iphs=1,nphases), &
					&	(yalpha(iphs),iphs=1,nphases), 			&
					&	norm_drag_chem, norm_drag_chem_spec(1:nphases)
				endif
				close (unitnum)
				close (unitnum_chem)

				filename1=trim(run_name)//"_drag_parts"
				unitnum = 1
				call instant_file_opener(filename1,unitnum,.true.)

				tmp1 = sqrt(dot_product(pres_total,pres_total)) / nbody
				tmp2 = sqrt(dot_product(visc_total,visc_total)) / nbody
				tmp3 = sqrt(dot_product(mpg,mpg)) * pi*dia_phys**3/6
				tmp4 = 3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*dia_phys

				re = (mixmeanslipmod+SMALL_NUMBER)*dia_phys / vis

				write (unitnum,"(5f10.4,4D15.7)") dbydx, lybyd, re_out, archno, mean_volfrac, &
				& 	tmp1/tmp4, tmp2/tmp4, tmp3/tmp4, (tmp1+tmp2+tmp3)/tmp4
				close(unitnum)
			else
				if (nphases==1) then
					nvar  = 9
					nvar1 = 5
					nvar2 = 4
				else
					nvar  = 12
					nvar1 = 9
					nvar2 = 3
				endif

				filename1 = "_drag_post.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line, include_nmis)

				if (nphases==1) then
					nvar  = 7
					nvar1 = 5
					nvar2 = 2
				else
					nvar  = 12
					nvar1 = 9
					nvar2 = 3
				endif
				filename1 = "_drag_chem_post.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line, include_nmis)

				nvar  = 9
				nvar1 = 5
				nvar2 = 4

				filename1 = "_drag_parts_post.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line, include_nmis)
			endif
		endif
	end subroutine write_drag


	SUBROUTINE compute_ibm_drag(phi, rem, F)
		IMPLICIT NONE
		Real(prcn), Intent(in) :: phi, rem
		Real(prcn), Intent(out) :: F
		Real(prcn) :: RE, FISOL, c, F0, F1

		RE =  Rem
		c = phi
		FISOL = 1.d0 + 0.15*(RE**0.687)
		FISOL = FISOL/(1-c)**3.d0
		F0 = 5.813*c/(1-c)**3.d0 + 0.485*c**(1.0/3.0)/(1-c)**4.d0
		F1 = RE*(c**3.d0)*(0.954 + 0.607*(c**3.d0/(1-c)**2.d0))
		F = FISOL + F0 + F1
	END SUBROUTINE compute_ibm_drag

	subroutine velocity_output
		implicit none
		character*100 filename1
		integer :: i, j, k, ii, jj, kk, idim, unitnum
		real(prcn) :: tmp, mean_energy, norm_vel
		

		unitnum = 1
		call outpufilename_gen(filename1,"velocity_field")
		open (unit=unitnum, file=trim(filename1), status="replace")
		write (unitnum,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
		write (unitnum,"(3(1A,1I),1a)") "i=", local_ni(1)+1, " j=", local_ni(2)+1, " k=", local_ni(3)+1, " f=point"

		if (zero_slip) then
			norm_vel = ucharmod
			if (iturbon) then
				mean_energy = tke_i
			elseif (Ret>small_number) then
				mean_energy = grant_i*3/2
			endif
		else
			norm_vel = sqrt(dot_product(ufmean(:)-usmean(:),ufmean(:)-usmean(:)))
			mean_energy = half*norm_vel**2
		endif

		do k=1, local_ni(3)+1
			do j=1, local_ni(2)+1
				do i=1, local_ni(1)+1
					ii = i + local_i_start(1)-1
					jj = j + local_i_start(2)-1
					kk = k + local_i_start(3)-1

					tmp = dot_product(ubcp(i,j,k,:),ubcp(i,j,k,:))/2
					write (unitnum,"(3i6,5d15.7)") ii, jj, kk, ubcp(i,j,k,:)/norm_vel ! tmp/mean_energy 
					!, ubcp(i,j,k,:)/norm_vel!, pbcp(i,j,k) / (half*mean_energy**2)
				enddo
			enddo
		enddo
		close (unitnum)

		call particle_snapshot
	end subroutine velocity_output

	subroutine particle_snapshot
		implicit none
		integer :: m, fileunit
		character*100 :: filename
		

		if (I_AM_NODE_ZERO) then
			fileunit=1
			filename = trim(run_name)//"_sphr_position.dat"
			open (unit=fileunit, file=trim(filename), status="replace")

			write (fileunit,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
			do m=1, nbody
				write(fileunit,'(10(2x,1d15.7))')  xc(m,:), radbdy(m)
			enddo
			close (fileunit)
		endif
	end subroutine particle_snapshot


	subroutine fluid_vel_cage
		use init_turb, only : mx_iso, u_prime
		implicit none
		integer :: index_in(ndim), index_out(ndim)
		real(prcn) :: position_in(ndim), position_out(ndim), l_cage
		integer :: i,j,k,ibody,imin,imax,jmin,jmax,kmin,kmax

		real(prcn), allocatable :: tmp_vel(:,:)
		integer, allocatable :: tmp_count(:), count_cage(:)

		if (I_AM_NODE_ZERO) write (*,*) "IN COMPUTE_UF_CAGE"
		if (allocated(uf_cage)) deallocate(uf_cage)
		allocate(uf_cage(nbody,ndim))

		allocate(tmp_vel(nbody,ndim))
		allocate(tmp_count(nbody))		
		allocate(count_cage(nbody))

		uf_cage = zero
		count_cage = 0
		do ibody=1, nbody
			!l_cage = dble(mx) / mx_iso * radbdy(ibody)*2
			l_cage = radbdy(ibody)*2 !* 2
			imin = floor( xc(ibody,1)-l_cage/2)
			jmin = floor( xc(ibody,2)-l_cage/2)
			kmin = floor( xc(ibody,3)-l_cage/2)

			imax = ceiling( xc(ibody,1)+l_cage/2)
			jmax = ceiling( xc(ibody,2)+l_cage/2)
			kmax = ceiling( xc(ibody,3)+l_cage/2)


			do k=kmin, kmax
				do j=jmin, jmax
					do i=imin, imax
						index_in(1) = i
						index_in(2) = j
						index_in(3) = k

						position_in(1) = i
						position_in(2) = j
						position_in(3) = k

						if (point_in_this_domain(index_in, index_out, position_in, position_out)) then
							if (fluid_atijk(index_out(1),index_out(2),index_out(3))) then
								count_cage(ibody) = count_cage(ibody) + 1
								uf_cage(ibody,:) = uf_cage(ibody,:) + ubcp(index_out(1),index_out(2),index_out(3),:)
							 endif
						endif
					enddo
				 enddo
			enddo
		enddo


		tmp_vel = uf_cage
		uf_cage = zero
		GLOBAL_DOUBLE_SUM(tmp_vel,uf_cage,nbody*ndim,comm_cart_2d)
		
		tmp_count = count_cage
		count_cage = 0
		GLOBAL_INT_SUM(tmp_count,count_cage,nbody,comm_cart_2d)

		do ibody=1, nbody
			uf_cage(ibody,:) = uf_cage(ibody,:)/count_cage(ibody)
		enddo

		deallocate(tmp_vel)
		deallocate(tmp_count)
		deallocate(count_cage)

		if (I_AM_NODE_ZERO) then
			write (*,*) "L_CAGE = ", l_cage, radbdy(1)
			write (*,*) "MAX UF_CAGE(1)", maxval(uf_cage(:,1)) / u_prime
			write (*,*) "MAX UF_CAGE(2)", maxval(uf_cage(:,2)) / u_prime
			write (*,*) "MAX UF_CAGE(3)", maxval(uf_cage(:,3)) / u_prime
		endif
	end subroutine fluid_vel_cage

	subroutine Apr_App_scatter
		use init_turb, only : u_prime
		implicit none

		real(prcn), allocatable :: particle_accl1(:,:), particle_accl2(:,:)
		real(prcn) :: tau_p, mass, tmp1, tmp2, norm
		integer :: ibody, unitnum
		character*100 filename

		call fluid_vel_cage


		tau_p = rhos/rhof / 18 * dia_phys**2/ vis
		if (I_AM_NODE_ZERO) write (*,*) "TAU_P = ", tau_p

		if (.not.allocated(particle_accl1)) allocate(particle_accl1(nbody,ndim))
		if (.not.allocated(particle_accl2)) allocate(particle_accl2(nbody,ndim))

		
		do ibody = 1, nbody
			particle_accl1(ibody,:) = (uf_cage(ibody,:)-velbdy(ibody,:)) / tau_p
		enddo

		do ibody = 1, nbody
			mass= rhos*pi*dia_phys**3.d0/6.d0
			particle_accl2(ibody,:) = force(ibody,:) / mass
		end do	

		norm = 18.d0 * rhof/rhos * vis / dia_phys * u_prime

		if (I_AM_NODE_ZERO) then

			call particle_reynolds_pdf


			unitnum = 1
			filename = trim(run_name)//"_Apr_App"
			call instant_file_opener(filename,unitnum,.true.)
			write (unitnum,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
			do ibody =1, nbody
				tmp1 = sqrt(dot_product(particle_accl1(ibody,:), particle_accl1(ibody,:)))
				tmp2 = sqrt(dot_product(particle_accl2(ibody,:), particle_accl2(ibody,:)))
				write (unitnum,'(500(2x, e20.12))') tmp1/norm, tmp2/norm
			enddo
			close(unitnum)
		endif
	end subroutine Apr_App_scatter


	subroutine compute_energy_spectrum_function
		use init_turb, only : statistics
		implicit none
		real(prcn) :: rmax(ndim), dist(ndim), length
		integer :: i, j, k, ii, jj, kk, idim, jdim, nparts, ibody, part_select, turn

		if (I_AM_NODE_ZERO) write (*,*) "IN compute_energy_spectrum_function"

		rmax(:) = dble(global_n(:))/2
		do turn=1, 1
			if (I_AM_NODE_ZERO) write (*,*) "PASS NUM = ", turn
			!if (turn == 2) ubcp(:,:,:,:) = zero

  			do idim=1, ndim
				urtmp(:,:,:) = zero
				do k=1, local_ni(3)
					do j=1, local_ni(2)
						do i=1, local_ni(1)
							if (.not.fluid_atijk(i,j,k)) then
								part_select = -1

								nparts = gnacc_part(i,j,k,1)
								if (nparts==0) then
									write (*,"(4i,1a)") myid, i, j, k, "A SOLID POINT NOT ASSOCIATED WITH ANY PARTICLE"
									PARALLEL_FINISH()
									stop
								elseif (nparts>1) then
									ii = i+local_i_start(1)-1
									jj = j+local_i_start(2)-1
									kk = k+local_i_start(3)-1
									do ibody=1, nparts
										dist(1) = abs(ii-xc( gnacc_part(i,j,k,ibody+1) ,1))
										dist(2) = abs(jj-xc( gnacc_part(i,j,k,ibody+1) ,2))
										dist(3) = abs(kk-xc( gnacc_part(i,j,k,ibody+1) ,3))
										do jdim=1, ndim
											if (dist(jdim)> rmax(jdim)) dist(jdim) = 2*rmax(jdim) - dist(jdim)
										enddo
										length = sqrt(dot_product(dist(:),dist(:)))
										if (length<=radbdy(gnacc_part(i,j,k,ibody+1))) then
											part_select = gnacc_part(i,j,k,ibody+1)
											exit
										endif
									enddo
								elseif (nparts==1) then
									part_select = gnacc_part(i,j,k,2)
								endif

								if (part_select==-1) then
									write (*,"(4i,1a)") myid, i, j, k, "NO PARTICLE WAS SELECTED FOR THIS SOLID POINT"
									PARALLEL_FINISH()
									stop
								endif

								if (turn==1) then
									urtmp(i,j,k) = velbdy(part_select,idim)
								else
									urtmp(i,j,k) = velbdy(part_select,idim)
								endif
							else
								if (turn==1) then
									urtmp(i,j,k) = ubcp(i,j,k,idim)
								else
									urtmp(i,j,k) = zero
								endif
							endif
						enddo
					enddo
				enddo

				call fftwr2c(urtmp, uftmp)

				if (turn==1) then
					u(:,:,:,idim) = uftmp(:,:,:)
				elseif (turn==2) then
					u(:,:,:,idim) = u(:,:,:,idim) - uftmp(:,:,:)
				endif
			enddo
		enddo

		call statistics(1)

	end subroutine compute_energy_spectrum_function


	subroutine outpufilename_gen(filename,var_name)
		implicit none

		Character(LEN=*),Intent(out) :: filename
		Character(LEN=*),Intent(in) :: var_name

		CHARACTER*10 :: filename1
		Character*100 :: filenameloc
		Integer :: node, strlen

		BARRIER(comm_cart_2d)

		if(I_AM_NODE_ZERO)then
#if PARALLEL
			!^^^^ MODIFIED FOR MORE PROCS
			if (nproc<=100) then
				write (filename1,fmt="('_NODE',1i2.2)") myid
			elseif (nproc<=1000) then
				write (filename1,fmt="('_NODE',1i3.3)") myid
			elseif (nproc<=10000) then
				write (filename1,fmt="('_NODE',1i4.4)") myid
			elseif (nproc<=100000) then
				write (filename1,fmt="('_NODE',1i5.5)") myid
			elseif (nproc<=1000000) then
				write (filename1,fmt="('_NODE',1i6.6)") myid
			endif
#else
			write (filename1,"(1a)") ""
#endif
			filename = TRIM(run_name)//'_'//TRIM(var_name)//TRIM(filename1)//'.dat'
			filenameloc = ""
		else
			filename = ""
		endif

		if(I_AM_NODE_ZERO)then
			do node=1,nproc-1
				!^^^^ MODIFIED FOR MORE PROCS
				if (nproc<=100) then
					write (filename1,fmt="('_NODE',1i2.2)") node
				elseif (nproc<=1000) then
					write (filename1,fmt="('_NODE',1i3.3)") node
				elseif (nproc<=10000) then
					write (filename1,fmt="('_NODE',1i4.4)") node
				elseif (nproc<=100000) then
					write (filename1,fmt="('_NODE',1i5.5)") node
				elseif (nproc<=1000000) then
					write (filename1,fmt="('_NODE',1i6.6)") node
				endif

				filenameloc = TRIM(run_name)//'_'//TRIM(var_name)//TRIM(filename1)//'.dat'
				strlen = LEN_trim(filenameloc)

				SEND_INT(strlen,1,node,0,comm_cart_2d)
				SEND_CHARACTER(filenameloc,strlen,node,1,comm_cart_2d)
			enddo
		else  
			!RECV_STRING(filename,strlen,node_zero,0,1,comm_cart_2d,status)
			RECV_INT(strlen,1,node_zero,0,comm_cart_2d,status)
			RECV_CHARACTER(filename,strlen,node_zero,1,comm_cart_2d,status)
		endif

		BARRIER(comm_cart_2d)
	end subroutine outpufilename_gen



	subroutine read_part_info(unitnum,count,io)
		implicit none
		integer, intent(in) :: unitnum
		integer, intent(inout) :: count, io
		character*80 filename
		logical :: filexist, isopen
		integer, save :: turn = 0
		real(prcn), save :: t_tmp

10		continue
		if (turn==12) goto 20

		if (turn==0) then
			filename = TRIM(RUN_NAME)//"_part_info.rst"
		elseif (turn==1) then
			filename = TRIM(RUN_NAME)//"_part_info_0.dat"
		elseif (turn==2) then
			filename = TRIM(RUN_NAME)//"_part_info_1.dat"
		elseif (turn==3) then
			filename = TRIM(RUN_NAME)//"_part_info_2.dat"
		elseif (turn==4) then
			filename = TRIM(RUN_NAME)//"_part_info_3.dat"
		elseif (turn==5) then
			filename = TRIM(RUN_NAME)//"_part_info_4.dat"
		elseif (turn==6) then
			filename = TRIM(RUN_NAME)//"_part_info_5.dat"
		elseif (turn==7) then
			filename = TRIM(RUN_NAME)//"_part_info_6.dat"
		elseif (turn==8) then
			filename = TRIM(RUN_NAME)//"_part_info_7.dat"
		elseif (turn==9) then
			filename = TRIM(RUN_NAME)//"_part_info_8.dat"
		elseif (turn==10) then
			filename = TRIM(RUN_NAME)//"_part_info_9.dat"
		elseif (turn==11) then
			filename = TRIM(RUN_NAME)//"_part_info.dat"
		endif

		inquire (file=trim(filename), exist=filexist,opened=isopen)
		if (filexist .and. .not.isopen) then
			write (*,*) "FILE "//filename//" EXISTS"
			open (unit=unitnum, file=trim(filename), status="old", action="read", form="unformatted")
			isopen = .true.
		elseif (.not.filexist) then
			write (*,*) "FILE  "//filename//" DOES NOT EXIST"
			turn = turn+1
			goto 10
		endif

		if (isopen) then
30			continue
			read(unitnum,iostat=io) t
			read(unitnum,iostat=io) xc(1:nbody, 1:ndim)
			read(unitnum,iostat=io) velbdy(1:nbody, 1:ndim)
			read(unitnum,iostat=io) force(1:nbody, 1:ndim)
			read(unitnum,iostat=io) pres(1:nbody, 1:ndim)
			read(unitnum,iostat=io) visc(1:nbody, 1:ndim)
			read(unitnum,iostat=io) contact_force(1:nbody, 1:ndim)
			read(unitnum,iostat=io) frame_vel(1:ndim)
			read(unitnum,iostat=io) frame_accln(1:ndim)
			read(unitnum,iostat=io) ufmean(1:ndim)

			if (io<0) then
				write (*,*) "END OF FILE ", trim(filename)
				close(unitnum)
				turn = turn+1
				goto 10
			else
				if (t<=t_tmp) goto 30
				t_tmp = t
				count = count+1
			endif
		endif
20		continue
	end subroutine read_part_info


	subroutine read_tke(unitnum,count,io)
		use init_turb, only : u_prime
		implicit none
		integer, intent(in) :: unitnum
		integer, intent(inout) :: count, io
		character*80 filename
		logical :: filexist, isopen
		real(prcn), allocatable :: tke_array(:)
		real(prcn), save :: t_tmp

		filename = TRIM(RUN_NAME)//"_tke.dat"
		if (.not.allocated(tke_array)) allocate(tke_array(17))

		inquire (file=trim(filename), exist=filexist,opened=isopen)
		if (filexist .and. .not.isopen) then
			write (*,*) "FILE "//filename//" EXISTS"
			open (unit=unitnum, file=trim(filename), status="old", action="read", form="formatted")
			isopen = .true.
		elseif (.not.filexist) then
			write (*,*) "FILE  "//filename//" DOES NOT EXIST"
		endif

		if (isopen) then
10			continue
			read(unitnum,*,iostat=io) tke_array(1:17)
			tke = tke_array(2)
			u_prime = sqrt(two*tke/3)

			if (io<0) then
				write (*,*) "END OF FILE ", trim(filename)
				close(unitnum)
			else
				if (t<=t_tmp) goto 10
				t_tmp = t
				count = count+1
			endif
		endif
	end subroutine read_tke


	subroutine kp_source_dissipation
		use init_turb, only : calc_initial_turb, eddy_time_i, tke_i, epsf_i
		implicit none
		integer :: i, count1, count2, io1, io2, partunit, tkeunit
		integer :: nvar, nvar1, nvar2
		character*100 filename


		if (iturbon) then
			call calc_initial_turb
			t_conv = eddy_time_i
		endif

		if (I_AM_NODE_ZERO) then
			write (*,*) "IN AiVi..."
			if (.not.post_no_flow_mem_alloc) then
				first_pass_post = .true.
				count1 = 0
				count2 = 0
				partunit = 10001
				!tkeunit  = 10002
				do
					call read_part_info(partunit,count1,io1)
					!call read_tke(tkeunit,count2,io2)
					!if (io1<0 .or. io2<0) exit

					if (io1<0) exit

					if (mod(count1,skip_num)==0) then
						write (*,*) "AT TIME = ", t/t_conv

						if (count1>1) then
							call compute_AIVI
							call particle_acceleration_pdf
							!call particle_reynolds_pdf
							first_pass_post = .false.
						endif
					endif
				enddo

				p_acc_avg(:) = p_acc_avg(:)/p_acc_avg_count

				filename = trim(run_name)//"_pdf_p_acc_avg.dat"
#if 0
				open (unit=1, file=trim(filename), status="replace", action="write")
				write (1,*) "zone"
				do i=1, nbins
					write (1,"(2d15.7)") radbin(i), p_acc_avg(i)
				enddo
				close(1)
#endif
			else
				nvar  = 8
				nvar1 = 6
				nvar2 = 2

				filename = "_AiVi_avg.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line, include_nmis)

				nvar  = 2
				nvar1 = 1
				nvar2 = 1

				line = nbins

				filename = "_pdf_p_acc.dat"
				call mis_average(nvar, nvar1, nvar2, filename, nbins, .false.)
			endif
		endif

	end subroutine kp_source_dissipation



	subroutine compute_AiVi
		implicit none
		real(prcn), allocatable :: acc_fluc(:,:), vel_fluc(:,:), acc_var(:,:), vel_var(:,:), acc_fluc_meanf(:,:), acc_var_meanf(:,:)

		allocate(acc_fluc(nbody,ndim), vel_fluc(nbody,ndim), acc_var(nphases,ndim), vel_var(nphases,ndim))
		allocate(acc_fluc_meanf(nbody,ndim), acc_var_meanf(nphases,ndim))

		call AiVi(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)

		deallocate(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)
	end subroutine compute_AiVi


	subroutine AiVi(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)
		use init_turb, only : eddy_time_i, tke_i, epsf_i, u_eta_i, tau_eta_i
		implicit none

		real(prcn), intent(out) :: acc_fluc(nbody,ndim), vel_fluc(nbody,ndim), acc_var(nphases,ndim), vel_var(nphases,ndim), acc_fluc_meanf(nbody,ndim), acc_var_meanf(nphases,ndim)

		integer :: pstart, pend, m, iphs, idim
		character*50 filename
		real(prcn), allocatable :: acc_avg(:,:), acc_avg_meanf(:,:), phase_mass(:)
		real(prcn) :: tmp1, re_tmp, f_tmp, norm_factor, norm_acc_factor, norm_vel_factor
		integer :: nvar1, nvar2, nvar
		logical :: filexist


		if (.not.allocated(acc_mean)) allocate(acc_mean(nphases), acc_variance(nphases), vel_mean(nphases))
		if (.not.allocated(acc_avg_meanf)) allocate(acc_avg_meanf(nphases,ndim))
		if (.not.allocated(acc_avg)) allocate(acc_avg(nphases,ndim), phase_mass(nphases))

#if 0
		! PARTICLE ACCELERATION AND ITS SD OBTAINED FROM THE MEAN DRAG FORCE MODEL
		acc_avg_meanf = 0
		pstart = 1
		do iphs = 1, nphases
			phase_mass(iphs) = rhos*pi*(phase_array(iphs)%dia)**3.d0/6.d0
			pend = pstart + phase_array(iphs)%npart- 1
			do m = pstart, pend
				re_tmp = sqrt(dot_product(ufmean(:)-velbdy(m,:), ufmean(:)-velbdy(m,:)))
				re_tmp = re_tmp * (1-mean_volfrac) * dia_phys / vis

				call compute_ibm_drag(mean_volfrac, re_tmp, f_tmp)
				acc_fluc_meanf(m,:) = f_tmp * 3.d0 * pi * vis * dia_phys * (1-mean_volfrac) * (ufmean(:)-velbdy(m,:)) / phase_mass(iphs)
				acc_avg_meanf(iphs,:) = acc_avg_meanf(iphs,:) + acc_fluc_meanf(m,:)
			enddo
			acc_avg_meanf(iphs,:) = acc_avg_meanf(iphs,:) / (pend-pstart+1)
			pstart = pend + 1
		end do	

		acc_var_meanf = zero

		do idim=1, ndim
			pstart = 1
			do iphs = 1, nphases
				phase_mass(iphs) = rhos*pi*(phase_array(iphs)%dia)**3.d0/6.d0
				pend = pstart + phase_array(iphs)%npart - 1
				do m = pstart, pend
					acc_fluc_meanf(m,idim) = acc_fluc_meanf(m,idim) - acc_avg_meanf(iphs,idim)

					acc_var_meanf(iphs,idim) = acc_var_meanf(iphs,idim) + acc_fluc_meanf(m,idim)**2
				enddo
				acc_var_meanf(iphs,idim) = acc_var_meanf(iphs,idim) / (pend-pstart+1)
				acc_var_meanf(iphs,idim) = sqrt(acc_var_meanf(iphs,idim))
				pstart = pend + 1
			enddo
		enddo
#endif


		if (zero_slip) then
			if (iturbon) then
				norm_acc_factor = u_eta_i/tau_eta_i
				norm_vel_factor = u_eta_i
			endif
		else
			norm_acc_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length) / phase_mass(iphs)
			norm_vel_factor = (mixmeanslipmod+SMALL_NUMBER)*maxvolfrac
		endif



		!MEAN ACCELERATION
		acc_avg = zero
		pstart = 1
		do iphs = 1, nphases
			phase_mass(iphs) = rhos*pi*(phase_array(iphs)%dia)**3.d0/6.d0
			pend = pstart + phase_array(iphs)%npart- 1
			do m = pstart, pend
				acc_avg(iphs,:) = acc_avg(iphs,:)+force(m,:)/phase_mass(iphs)/norm_acc_factor
			enddo
			acc_avg(iphs,:) = acc_avg(iphs,:) / (pend-pstart+1)
			pstart = pend + 1
		enddo

		do iphs = 1, nphases
			acc_avg(iphs,:) = acc_avg(iphs,:) ! + mpg(:)*pi*((two*radbdy(1)*dx)**3.d0)/6.d0/phase_mass(iphs)
		enddo

		do iphs = 1, nphases
			acc_mean(iphs) = sqrt(dot_product(acc_avg(iphs,:), acc_avg(iphs,:)))
			vel_mean(iphs) = sqrt(dot_product(phase_array(iphs)%mean_spec_vel(:), phase_array(iphs)%mean_spec_vel(:)))
		enddo

		!FLUCTUATING ACCELERATION
		acc_fluc = zero
		pstart = 1
		do iphs = 1, nphases
			pend = pstart + phase_array(iphs)%npart- 1
			do m = pstart, pend
				do idim=1, ndim
					acc_fluc(m,idim) = force(m,idim)/phase_mass(iphs)/norm_acc_factor - acc_avg(iphs,idim)
				end do
			enddo
			pstart = pend + 1
		end do

		!FLUCTUATING VELOCITY
		vel_fluc = zero
		pstart = 1
		do iphs = 1, nphases
			pend = pstart + phase_array(iphs)%npart- 1
			do m = pstart, pend
				do idim=1, ndim
					vel_fluc(m,idim) = (velbdy(m,idim) - phase_array(iphs)%mean_spec_vel(idim)) / norm_vel_factor
				end do
			enddo
			pstart = pend + 1
		end do

		!STANDARD DEVIATIONS
		acc_var = zero
		vel_var = zero

		do idim=1, ndim
			pstart = 1
			do iphs = 1, nphases
				pend = pstart + phase_array(iphs)%npart- 1
				do m = pstart, pend
					acc_var(iphs,idim) = acc_var(iphs,idim) + (force(m,idim)/phase_mass(iphs)/norm_acc_factor - acc_avg(iphs,idim))**2.d0
					vel_var(iphs,idim) = vel_var(iphs,idim) + ((velbdy(m,idim) - phase_array(iphs)%mean_spec_vel(idim)) / norm_vel_factor)**2.d0
				enddo

				acc_var(iphs,idim) = acc_var(iphs,idim) / (pend-pstart+1) !/ ndim
				vel_var(iphs,idim) = vel_var(iphs,idim) / (pend-pstart+1) !/ ndim

				acc_var(iphs,idim) = sqrt(acc_var(iphs,idim))
				vel_var(iphs,idim) = sqrt(vel_var(iphs,idim))

				pstart = pend + 1
			end do
		enddo

		do iphs=1, nphases
			acc_variance(iphs) = sqrt(dot_product(acc_var(iphs,:),acc_var(iphs,:)))
		enddo


		!SOURCE / DISSIPATION
		source = zero
		dissipation = zero
		pstart = 1
		do iphs = 1, nphases
			pend = pstart + phase_array(iphs)%npart- 1
			do m = pstart, pend
				tmp1 = dot_product(acc_fluc(m,:), vel_fluc(m,:))

				if (tmp1>zero) then
					source      = source      + tmp1
				else
					dissipation = dissipation - tmp1
				endif
			enddo
			pstart = pend + 1
		end do
		source      = source      / real(nbody, prcn) !*two/three 
		dissipation = dissipation / real(nbody, prcn) !*two/three 


		filename=trim(run_name)//"_AiVi.dat"

		if (first_pass_post) then	
			open (unit=1, file=trim(filename), status="replace", action="write")
		else
			inquire (file=trim(filename), exist=filexist)
			if (.not.filexist) then
				open (unit=1, file=trim(filename), status="replace", action="write")
			else
				open (unit=1, file=trim(filename), status="old", action="write", position="append")
			endif
		endif

		write (1,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
		do idim=1, ndim
			pstart = 1
			do iphs = 1, nphases
				pend = pstart + phase_array(iphs)%npart- 1
				do m = pstart, pend

					!if (vel_var(iphs,idim)>post_small.and.acc_var(iphs,idim)>post_small) then
					write (1,"(3D15.7)") acc_fluc(m,idim)/acc_var(iphs,idim), &
						&						vel_fluc(m,idim)/vel_var(iphs,idim) !, &
								!&	acc_fluc(m,idim)*vel_fluc(m,idim) / acc_var(iphs,idim)/vel_var(iphs,idim)
					!endif
				enddo
				pstart = pend + 1
			enddo
		enddo

!				do idim=1, ndim
!					write (1,*) "zone"
!					pstart = 1
!					do iphs = 1, nphases
!						pend = pstart + phase_array(iphs)%npart- 1
!						do m = pstart, pend
!!							if (acc_var_meanf(iphs,idim)>post_small.and.vel_var(iphs,idim)>post_small) then
!								write (1,"(3D15.7)") acc_fluc_meanf(m,idim)/acc_var_meanf(iphs,idim), vel_fluc(m,idim)/vel_var(iphs,idim), &
!									&	acc_fluc_meanf(m,idim)*vel_fluc(m,idim) / acc_var_meanf(iphs,idim)/vel_var(iphs,idim)
!!							endif
!						enddo
!						pstart = pend + 1
!					end do
!				enddo
		close (1)


		filename = trim(run_name)//"_AiVi_sd.dat"
		if (first_pass_post) then
			open (unit=2, file=trim(filename), status="replace", action="write")
			write (2,*) "zone"
		else
			inquire (file=trim(filename), exist=filexist)
			if (.not.filexist) then
				open (unit=2, file=trim(filename), status="replace", action="write")
				write (2,*) "zone"
			else
				open (unit=2, file=trim(filename), status="old", action="write", position="append")
			endif
		endif

		write (2,"(7D15.7)") t/t_conv, acc_var(1,:), vel_var(1,:)
		!endif
		close (2)


		filename = trim(run_name)//"_AiVi_avg.dat"
		if (first_pass_post) then
			open (unit=1, file=trim(filename), status="replace", action="write")
			write (1,*) "zone"
		else
			inquire (file=trim(filename), exist=filexist)
			if (.not.filexist) then
				open (unit=1, file=trim(filename), status="replace", action="write")
				write (1,*) "zone"
			else
				open (unit=1, file=trim(filename), status="old", action="write", position="append")
			endif
		endif

		tmp1 = sqrt( dot_product (acc_var(1,:), vel_var(1,:)) )

		if (zero_slip) then
			if (iturbon) then
				norm_factor = epsf_i
			endif
		endif

		!if (from_post) then
		!	write (1,"(6f10.3,2D15.7)") dbydx, lybyd, mean_volfrac, re, rhos/rhof, coeff_rest, source/tmp1, dissipation/tmp1
		!else
			!write (1,"(3D15.7)") t/t_conv, source/tmp1, dissipation/tmp1
			write (1,"(4D15.7)") t/t_conv, source/norm_factor, dissipation/norm_factor, (source-dissipation)/norm_factor
		!endif
		close (1)
		write (*,"(1A,2d15.7)") "SOURCE, DISSIPATION = ", source/tmp1, dissipation/tmp1
	end subroutine AiVi


	subroutine particle_acceleration_pdf
		use init_turb, only : u_eta_i, tau_eta_i
		use bcsetarrays, only : ppr, diffn
		implicit none

		real(prcn), allocatable :: p_acc(:,:), hist(:), radbin(:), tmp_array(:)
		real(prcn) :: u_min, u_max
		integer :: ibody, idim, i, j, k, count, unit1, nvar, nvar1, nvar2, ip, im
		character*50 filename1, filename2, filename3, filename4, filename6, filename7, filename8, filename9
		real(prcn) :: mean1, mean2, mean3, mean4, mean_tmp
		real(prcn) :: sd1, sd2, sd3, sd4, sd_tmp, norm_factor, mass
		logical :: filexist, in_parallel


		mass = rhos * pi/6 * dia_phys**3
		if (zero_slip) then
			if (iturbon) then
				norm_factor = mass * u_eta_i/tau_eta_i
			endif
		else
			norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)
		endif

		if (.not.allocated(p_acc_avg)) then
			allocate(p_acc_avg(nbins))
			p_acc_avg = zero
			p_acc_avg_count = 0
		endif
    



		filename1 = trim(run_name)//"_pdf_p_acc1.dat"
		filename2 = trim(run_name)//"_pdf_p_acc2.dat"
		filename3 = trim(run_name)//"_pdf_p_acc3.dat"
		filename4 = trim(run_name)//"_pdf_p_accT.dat"

		filename6 = trim(run_name)//"_pdf_p_mean_sd1.dat"
		filename7 = trim(run_name)//"_pdf_p_mean_sd2.dat"
		filename8 = trim(run_name)//"_pdf_p_mean_sd3.dat"
		filename9 = trim(run_name)//"_pdf_p_mean_sdT.dat"


		if (first_pass_post) then
			open (unit=1, file=trim(filename1), status="replace", action="write")
			open (unit=2, file=trim(filename2), status="replace", action="write")
			open (unit=3, file=trim(filename3), status="replace", action="write")
			open (unit=4, file=trim(filename4), status="replace", action="write")

			open (unit=6, file=trim(filename6), status="replace", action="write")
			open (unit=7, file=trim(filename7), status="replace", action="write")
			open (unit=8, file=trim(filename8), status="replace", action="write")
			open (unit=9, file=trim(filename9), status="replace", action="write")
		else
			inquire (file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				open (unit=1, file=trim(filename1), status="replace", action="write")
			else
				open (unit=1, file=trim(filename1), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename2), exist=filexist)
			if (.not.filexist) then
				open (unit=2, file=trim(filename2), status="replace", action="write")
			else
				open (unit=2, file=trim(filename2), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename3), exist=filexist)
			if (.not.filexist) then
				open (unit=3, file=trim(filename3), status="replace", action="write")
			else
				open (unit=3, file=trim(filename3), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename4), exist=filexist)
			if (.not.filexist) then
				open (unit=4, file=trim(filename4), status="replace", action="write")
			else
				open (unit=4, file=trim(filename4), status="old", action="write", position="append")
			endif


			inquire (file=trim(filename6), exist=filexist)
			if (.not.filexist) then
				open (unit=6, file=trim(filename6), status="replace", action="write")
			else
				open (unit=6, file=trim(filename6), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename7), exist=filexist)
			if (.not.filexist) then
				open (unit=7, file=trim(filename7), status="replace", action="write")
			else
				open (unit=7, file=trim(filename7), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename8), exist=filexist)
			if (.not.filexist) then
				open (unit=8, file=trim(filename8), status="replace", action="write")
			else
				open (unit=8, file=trim(filename8), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename9), exist=filexist)
			if (.not.filexist) then
				open (unit=9, file=trim(filename9), status="replace", action="write")
			else
				open (unit=9, file=trim(filename9), status="old", action="write", position="append")
			endif
		endif

		write (1,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
		write (2,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
		write (3,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
		write (4,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '

		allocate(p_acc(nbody, ndim))
		p_acc = zero

		do ibody=1, nbody
			p_acc(ibody,:) = force(ibody,:)  / norm_factor
		enddo

		allocate(hist(nbins), radbin(nbins))
		allocate(tmp_array(nbody))

		in_parallel = .false.
		do idim=1, ndim
			tmp_array(:) = p_acc(:,idim)
			call make_histogram(tmp_array, nbody, nbody, nbins, radbin, hist, mean_tmp, sd_tmp, pdf_normalized , in_parallel)
			if (idim==1) then
				mean1 = mean_tmp
				sd1 = sd_tmp
			elseif (idim==2) then
				mean2 = mean_tmp
				sd2 = sd_tmp
			elseif (idim==3) then
				mean3 = mean_tmp
				sd3 = sd_tmp
			endif

			do i=1, nbins
				!if (hist(i)>small_number)
				write (idim,"(2d15.7)") radbin(i), hist(i)
			enddo

			if (idim==1) then
				p_acc_avg(:) = hist(:)
				p_acc_avg_count = p_acc_avg_count+1
			endif

		enddo

		do i=1, nbody
			tmp_array(i) = sqrt(dot_product(p_acc(i,:),p_acc(i,:)))
		enddo

		call make_histogram(tmp_array, nbody, nbody, nbins, radbin, hist, mean_tmp, sd_tmp, pdf_normalized, in_parallel)

		mean4 = mean_tmp
		sd4 = sd_tmp

		do i=1, nbins
			!if (hist(i)>small_number)
			write (4,"(2d15.7)") radbin(i), hist(i)
		enddo

		write (6,"(3d15.7)") t/t_conv, mean1, sd1
		write (7,"(3d15.7)") t/t_conv, mean2, sd2
		write (8,"(3d15.7)") t/t_conv, mean3, sd3
		write (9,"(3d15.7)") t/t_conv, mean4, sd4

		close (1)
		close (2)
		close (3)
		close (4)

		close (6)
		close (7)
		close (8)
		close (9)
	end subroutine particle_acceleration_pdf


	subroutine particle_reynolds_pdf
		use init_turb, only : u_prime
		implicit none

		real(prcn), allocatable :: p_vel(:,:), hist(:), radbin(:), tmp_array(:)
		real(prcn) :: u_min, u_max
		integer :: ibody, idim, i, j, k, count, unit1, nvar, nvar1, nvar2, ip, im
		character*50 filename1, filename2, filename3, filename4, filename6, filename7, filename8, filename9
		real(prcn) :: mean1, mean2, mean3, mean4, mean_tmp
		real(prcn) :: sd1, sd2, sd3, sd4, sd_tmp, norm_factor
		logical :: filexist, in_parallel


		if (zero_slip) then
			if (iturbon) then
				norm_factor = u_prime
			endif
		else
			norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)
		endif



		filename1 = trim(run_name)//"_pdf_Rep.dat"

		filename2 = trim(run_name)//"_pdf_Rep_mean_sd.dat"


		if (first_pass_post) then
			open (unit=1, file=trim(filename1), status="replace", action="write")
			open (unit=2, file=trim(filename2), status="replace", action="write")
		else
			inquire (file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				open (unit=1, file=trim(filename1), status="replace", action="write")
			else
				open (unit=1, file=trim(filename1), status="old", action="write", position="append")
			endif

			inquire (file=trim(filename2), exist=filexist)
			if (.not.filexist) then
				open (unit=2, file=trim(filename2), status="replace", action="write")
			else
				open (unit=2, file=trim(filename2), status="old", action="write", position="append")
			endif
		endif

		write (1,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '

		allocate(p_vel(nbody, ndim))
		p_vel = zero

		do ibody=1, nbody
			p_vel(ibody,:) = velbdy(ibody,:)
		enddo

		allocate(hist(nbins), radbin(nbins))
		allocate(tmp_array(nbody))

!		do ibody=1, nbody
!			tmp_array(ibody) = abs(sqrt(dot_product(p_vel(ibody,:),p_vel(ibody,:)))-u_prime) * dia_phys / vis
!		enddo

		do ibody=1, nbody
			tmp_array(ibody) = sqrt(dot_product(uf_cage(ibody,:)-velbdy(ibody,:),uf_cage(ibody,:)-velbdy(ibody,:))) * dia_phys / vis
		enddo


		in_parallel = .false.
		call make_histogram(tmp_array, nbody, nbody, nbins, radbin, hist, mean1, sd1, .false., in_parallel)

		write (2,"(3d15.7)") t/t_conv, mean1, sd1

		close (1)
		close (2)
	end subroutine particle_reynolds_pdf




	subroutine post_part_stat
		use init_turb, only : epsf_forcing, forced_turbulence, eddy_time_i, calc_initial_turb
		implicit none
		integer :: ibody, count1, io1, partunit
		integer :: nvar, nvar1, nvar2
		character*100 filename


		if (iturbon) then
			call calc_initial_turb
			t_conv = eddy_time_i
		endif

		if (I_AM_NODE_ZERO) then
			write (*,*) "IN POST_PART_STAT..."
			if (.not.post_no_flow_mem_alloc) then
				first_pass_post = .true.
				count1 = 0
				partunit = 10001

				filename=trim(run_name)//"_part_position.dat"
				open (unit=1, file=trim(filename), status="replace", action="write")

				do
					call read_part_info(partunit,count1,io1)
					if (io1<0) exit

					if (mod(count1,skip_num)==0) then
						write (*,*) "AT TIME = ", t
						if (count1>1) then

							write (1,"(1a,1d15.7,1a)") 'ZONE T= "', t, ' " '  ! t/t_conv, ' " '
							do ibody=1, nbody
								write (1,"(4d15.7)") xc(ibody,:), radbdy(ibody)
							enddo
						endif
					endif
				enddo
				close(1)
			else


			endif
		endif

	end subroutine post_part_stat




	subroutine fluid_acceleration_pdf
		use init_turb, only : u_prime
		use bcsetarrays, only : ppr, diffn
		implicit none

		real(prcn), allocatable :: hist(:), radbin(:), tmp_array(:)
		real(prcn) :: u_min, u_max
		integer :: ibody, idim, i, j, k, count, unit1, nvar, nvar1, nvar2
		character*50 filename1, filename2, filename3, filename4, filename6, filename7, filename8, filename9
		real(prcn) :: mean1, mean2, mean3, mean4, mean_tmp
		real(prcn) :: sd1, sd2, sd3, sd4, sd_tmp, norm_factor
		logical :: filexist, in_parallel
		integer :: nsample, ntotal

		if (zero_slip) then
			if (iturbon) then
				norm_factor = u_prime
			endif
		  else
			norm_factor = mixmeanslipmod/(one-maxvolfrac)
		endif

		first_pass_post = .true.
		if (I_AM_NODE_ZERO) then
			write (*,*) "IN FLUID_ACCELERATION_PDF ..."
			filename1 = trim(run_name)//"_pdf_f_acc1.dat"
			filename2 = trim(run_name)//"_pdf_f_acc2.dat"
			filename3 = trim(run_name)//"_pdf_f_acc3.dat"
			filename4 = trim(run_name)//"_pdf_f_accT.dat"

			filename6 = trim(run_name)//"_pdf_f_mean_sd1.dat"
			filename7 = trim(run_name)//"_pdf_f_mean_sd2.dat"
			filename8 = trim(run_name)//"_pdf_f_mean_sd3.dat"
			filename9 = trim(run_name)//"_pdf_f_mean_sdT.dat"


			if (first_pass_post) then
				open (unit=1, file=trim(filename1), status="replace", action="write")
				open (unit=2, file=trim(filename2), status="replace", action="write")
				open (unit=3, file=trim(filename3), status="replace", action="write")
				open (unit=4, file=trim(filename4), status="replace", action="write")

				open (unit=6, file=trim(filename6), status="replace", action="write")
				open (unit=7, file=trim(filename7), status="replace", action="write")
				open (unit=8, file=trim(filename8), status="replace", action="write")
				open (unit=9, file=trim(filename9), status="replace", action="write")
			else
				inquire (file=trim(filename1), exist=filexist)
				if (.not.filexist) then
					open (unit=1, file=trim(filename1), status="replace", action="write")
				else
					open (unit=1, file=trim(filename1), status="old", action="write", position="append")
				endif

				inquire (file=trim(filename2), exist=filexist)
				if (.not.filexist) then
					open (unit=2, file=trim(filename2), status="replace", action="write")
				else
					open (unit=2, file=trim(filename2), status="old", action="write", position="append")
				endif

				inquire (file=trim(filename3), exist=filexist)
				if (.not.filexist) then
					open (unit=3, file=trim(filename3), status="replace", action="write")
				else
					open (unit=3, file=trim(filename3), status="old", action="write", position="append")
				endif

				inquire (file=trim(filename4), exist=filexist)
				if (.not.filexist) then
					open (unit=4, file=trim(filename4), status="replace", action="write")
				else
					open (unit=4, file=trim(filename4), status="old", action="write", position="append")
				endif


				inquire (file=trim(filename6), exist=filexist)
				if (.not.filexist) then
					open (unit=6, file=trim(filename6), status="replace", action="write")
				else
					open (unit=6, file=trim(filename6), status="old", action="write", position="append")
				endif

				inquire (file=trim(filename7), exist=filexist)
				if (.not.filexist) then
					open (unit=7, file=trim(filename7), status="replace", action="write")
				else
					open (unit=7, file=trim(filename7), status="old", action="write", position="append")
				endif

				inquire (file=trim(filename8), exist=filexist)
				if (.not.filexist) then
					open (unit=8, file=trim(filename8), status="replace", action="write")
				else
					open (unit=8, file=trim(filename8), status="old", action="write", position="append")
				endif

				inquire (file=trim(filename9), exist=filexist)
				if (.not.filexist) then
					open (unit=9, file=trim(filename9), status="replace", action="write")
				else
					open (unit=9, file=trim(filename9), status="old", action="write", position="append")
				endif
			endif

			write (1,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
			write (2,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
			write (3,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
			write (4,"(1a,1d15.7,1a)") 'ZONE T= "', t/t_conv, ' " '
		endif

		nsample = local_ni(1)*local_ni(2)*local_ni(3)
		ntotal  = global_n(1)*global_n(2)*global_n(3)

		allocate(hist(nbins), radbin(nbins))
		allocate(tmp_array(nsample))

		in_parallel = .true.
		urtmp(:,:,:) = zero
		do idim=1, ndim
			count = 0
			do k=1, local_ni(3)
				do j=1, local_ni(2)
					do i=1, local_ni(1)
						count = count+1
						tmp_array(count) = ubcp(i,j,k,idim) - ufmean(idim)
						urtmp(i,j,k) = urtmp(i,j,k) + (ubcp(i,j,k,idim) - ufmean(idim))**2
					enddo
				enddo
			enddo


			call make_histogram(tmp_array, nsample, ntotal, nbins, radbin, hist, mean_tmp, sd_tmp, pdf_normalized, in_parallel)
			if (idim==1) then
				mean1 = mean_tmp
				sd1 = sd_tmp
			elseif (idim==2) then
				mean2 = mean_tmp
				sd2 = sd_tmp
			elseif (idim==3) then
				mean3 = mean_tmp
				sd3 = sd_tmp
			endif

			if (I_AM_NODE_ZERO) then
				do i=1, nbins
					!if (hist(i)>small_number)
					write (idim,"(2d15.7)") radbin(i), hist(i)
				enddo
			endif
		enddo

		count = 0
		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					count = count+1
					tmp_array(count) = sqrt(urtmp(i,j,k))
				enddo
			enddo
		enddo

		call make_histogram(tmp_array, nsample, ntotal, nbins, radbin, hist, mean_tmp, sd_tmp, pdf_normalized, in_parallel)

		mean4 = mean_tmp
		sd4 = sd_tmp

		if (I_AM_NODE_ZERO) then
			do i=1, nbins
				!if (hist(i)>small_number)
				write (4,"(2d15.7)") radbin(i), hist(i)
			enddo


			write (6,"(3d15.7)") t/t_conv, mean1/norm_factor, sd1/norm_factor
			write (7,"(3d15.7)") t/t_conv, mean2/norm_factor, sd2/norm_factor
			write (8,"(3d15.7)") t/t_conv, mean3/norm_factor, sd3/norm_factor
			write (9,"(3d15.7)") t/t_conv, mean4/norm_factor, sd4/norm_factor

			close (1)
			close (2)
			close (3)
			close (4)

			close (6)
			close (7)
			close (8)
			close (9)
		endif
	end subroutine fluid_acceleration_pdf


	subroutine make_histogram(array, nsample, ntotal, nbins, radbin, hist, mean, sd, normalized, parallel_mode)
		implicit none
		integer, intent(in) :: nsample, ntotal, nbins
		real(prcn), intent(inout) :: array(nsample)
		real(prcn), intent(out) :: radbin(nbins), hist(nbins)
		real(prcn), intent(out) :: mean, sd
		logical, intent(in) :: normalized, parallel_mode

		real(prcn) :: var, left, right, dr, tmp
		real(prcn), allocatable :: hist_tmp(:)
		integer :: i, ibin
		real(prcn) :: sum1, sum2, sum3

		call calc_avr_var(nsample, ntotal, array, mean, var, parallel_mode)

		sd = sqrt(var/ntotal)

		if (normalized) array(:) = (array(:) - mean) / sd

		left = minval(array(:))
		right = maxval(array(:))

		if (parallel_mode) then
			tmp = left
			GLOBAL_DOUBLE_MIN(tmp,left,1,comm_cart_2d)

			tmp = right
			GLOBAL_DOUBLE_MAX(tmp,right,1,comm_cart_2d)
		endif

		dr = (right-left) / nbins

		do i=1, nbins
			radbin(i) = left + (i-.5)*dr
		enddo

		hist = zero
		do i=1, nsample
			ibin = (array(i)-left) / dr + 1
			if (ibin>nbins) ibin = nbins
			hist(ibin) = hist(ibin) + 1
		enddo

		if (parallel_mode) then
			allocate(hist_tmp(nbins))
			hist_tmp(:) = hist(:)
			hist(:) = zero_slip

			GLOBAL_DOUBLE_SUM(hist_tmp,hist,nbins,comm_cart_2d)
			deallocate(hist_tmp)
		endif

#if 0
		sum1=zero
		sum2=zero
		sum3=zero

		do i=1, nbins
			if (radbin(i)<=1) sum1 = sum1 + hist(i)
			if (1<radbin(i) .and. radbin(i)<=2) sum2 = sum2 + hist(i)
			if (2<radbin(i) .and. radbin(i)<=3) sum3 = sum3 + hist(i)
		enddo

		write (*,*) "PART SUM =", sum1/ ntotal
		write (*,*) "PART SUM =", sum2/ ntotal
		write (*,*) "PART SUM =", sum3/ ntotal
#endif

		hist(:) = hist(:) / ntotal / dr

		do i=1, nbins
			!if (hist(i)>small_number)
			write (1,"(2d15.7)") radbin(i), hist(i)
		enddo

	end subroutine make_histogram


	subroutine mis_average(n, n1, n2, string, nline, include_mis)
		implicit none

		integer, intent(in) :: n, n1, n2
		integer, intent(inout) :: nline
		character*50, intent(in) :: string
		logical, intent(in) :: include_mis

		integer :: imis,ivar, size, iline, io, nmis_tmp, nline_max=100000, nline_array(nmis)
		real(prcn), allocatable, dimension(:,:,:) :: var
		real(prcn), allocatable, dimension(:) :: avr_var, var_var
		real(prcn) :: confint
		character*50 filename, outform, outform_drag, zone
		character*2 tmp1, tmp2, tmp3, s1, s2, s3
		character*3 tmp_str, str
		logical :: filexist
		real(prcn) :: fisol, ff, ffre, ft

		nline_array = 0

		if (include_mis) then
			call to_string(nmis,str,size)
			tmp_str = ""
			do ivar=1, 3-size
				tmp_str="0"//trim(tmp_str)
			enddo
			tmp_str= trim(tmp_str)//trim(str)
		endif


		call to_string(n1,s1,size)

		tmp1 = ""
		do ivar=1, 2-size
			tmp1="0"//trim(tmp1)
		enddo
		tmp1=trim(tmp1)//trim(s1)
!		write (*,*) tmp1

		call to_string(n2,s2,size)

		tmp2 = ""
		do ivar=1, 2-size
			tmp2="0"//trim(tmp2)
		enddo
		tmp2=trim(tmp2)//trim(s2)
!		write (*,*) tmp2

		outform = "("//trim(tmp1)
		outform = trim(outform)//"f10.4,"
		outform = trim(outform)//trim(tmp2)
		outform = trim(outform)//"D15.7)"
		if (.not.include_mis) then
			outform = "(2"//trim(outform)//")"
		else
			outform = "(1i,2"//trim(outform)//")"
		endif

!		outform = trim(outform)//"D15.7)"
!		write (*,*) trim(outform)
!		read(*,*)
!		outform =       trim(outform)//")"


		nmis_tmp = 0
		do imis=1, nmis
			if (imis==1) then
				filename="MIS1"//trim(string)
			elseif (imis==2) then
				filename="MIS2"//trim(string)
			elseif (imis==3) then
				filename="MIS3"//trim(string)
			elseif (imis==4) then
				filename="MIS4"//trim(string)
			elseif (imis==5) then
				filename="MIS5"//trim(string)
			elseif (imis==6) then
				filename="MIS6"//trim(string)
			elseif (imis==7) then
				filename="MIS7"//trim(string)
			elseif (imis==8) then
				filename="MIS8"//trim(string)
			elseif (imis==9) then
				filename="MIS9"//trim(string)
			elseif (imis==10) then
				filename="MIS10"//trim(string)

			elseif (imis==11) then
				filename="MIS11"//trim(string)
			elseif (imis==12) then
				filename="MIS12"//trim(string)
			elseif (imis==13) then
				filename="MIS13"//trim(string)
			elseif (imis==14) then
				filename="MIS14"//trim(string)
			elseif (imis==15) then
				filename="MIS15"//trim(string)
			elseif (imis==16) then
				filename="MIS16"//trim(string)
			elseif (imis==17) then
				filename="MIS17"//trim(string)
			elseif (imis==18) then
				filename="MIS18"//trim(string)
			elseif (imis==19) then
				filename="MIS19"//trim(string)
			elseif (imis==20) then
				filename="MIS20"//trim(string)

			elseif (imis==21) then
				filename="MIS21"//trim(string)
			elseif (imis==22) then
				filename="MIS22"//trim(string)
			elseif (imis==23) then
				filename="MIS23"//trim(string)
			elseif (imis==24) then
				filename="MIS24"//trim(string)
			elseif (imis==25) then
				filename="MIS25"//trim(string)
			elseif (imis==26) then
				filename="MIS26"//trim(string)
			elseif (imis==27) then
				filename="MIS27"//trim(string)
			elseif (imis==28) then
				filename="MIS28"//trim(string)
			elseif (imis==29) then
				filename="MIS29"//trim(string)
			elseif (imis==30) then
				filename="MIS30"//trim(string)

			elseif (imis==31) then
				filename="MIS31"//trim(string)
			elseif (imis==32) then
				filename="MIS32"//trim(string)
			elseif (imis==33) then
				filename="MIS33"//trim(string)
			elseif (imis==34) then
				filename="MIS34"//trim(string)
			elseif (imis==35) then
				filename="MIS35"//trim(string)
			elseif (imis==36) then
				filename="MIS36"//trim(string)
			elseif (imis==37) then
				filename="MIS37"//trim(string)
			elseif (imis==38) then
				filename="MIS38"//trim(string)
			elseif (imis==39) then
				filename="MIS39"//trim(string)
			elseif (imis==40) then
				filename="MIS40"//trim(string)

			elseif (imis==41) then
				filename="MIS41"//trim(string)
			elseif (imis==42) then
				filename="MIS42"//trim(string)
			elseif (imis==43) then
				filename="MIS43"//trim(string)
			elseif (imis==44) then
				filename="MIS44"//trim(string)
			elseif (imis==45) then
				filename="MIS45"//trim(string)
			elseif (imis==46) then
				filename="MIS46"//trim(string)
			elseif (imis==47) then
				filename="MIS47"//trim(string)
			elseif (imis==48) then
				filename="MIS48"//trim(string)
			elseif (imis==49) then
				filename="MIS49"//trim(string)
			elseif (imis==50) then
				filename="MIS50"//trim(string)

			elseif (imis==51) then
				filename="MIS51"//trim(string)
			elseif (imis==52) then
				filename="MIS52"//trim(string)
			elseif (imis==53) then
				filename="MIS53"//trim(string)
			elseif (imis==54) then
				filename="MIS54"//trim(string)
			elseif (imis==55) then
				filename="MIS55"//trim(string)
			elseif (imis==56) then
				filename="MIS56"//trim(string)
			elseif (imis==57) then
				filename="MIS57"//trim(string)
			elseif (imis==58) then
				filename="MIS58"//trim(string)
			elseif (imis==59) then
				filename="MIS59"//trim(string)
			elseif (imis==60) then
				filename="MIS60"//trim(string)

			elseif (imis==61) then
				filename="MIS61"//trim(string)
			elseif (imis==62) then
				filename="MIS62"//trim(string)
			elseif (imis==63) then
				filename="MIS63"//trim(string)
			elseif (imis==64) then
				filename="MIS64"//trim(string)
			elseif (imis==65) then
				filename="MIS65"//trim(string)
			elseif (imis==66) then
				filename="MIS66"//trim(string)
			elseif (imis==67) then
				filename="MIS67"//trim(string)
			elseif (imis==68) then
				filename="MIS68"//trim(string)
			elseif (imis==69) then
				filename="MIS69"//trim(string)
			elseif (imis==70) then
				filename="MIS70"//trim(string)


			elseif (imis==71) then
				filename="MIS71"//trim(string)
			elseif (imis==72) then
				filename="MIS72"//trim(string)
			elseif (imis==73) then
				filename="MIS73"//trim(string)
			elseif (imis==74) then
				filename="MIS74"//trim(string)
			elseif (imis==75) then
				filename="MIS75"//trim(string)
			elseif (imis==76) then
				filename="MIS76"//trim(string)
			elseif (imis==77) then
				filename="MIS77"//trim(string)
			elseif (imis==78) then
				filename="MIS78"//trim(string)
			elseif (imis==79) then
				filename="MIS79"//trim(string)
			elseif (imis==80) then
				filename="MIS80"//trim(string)

			elseif (imis==81) then
				filename="MIS81"//trim(string)
			elseif (imis==82) then
				filename="MIS82"//trim(string)
			elseif (imis==83) then
				filename="MIS83"//trim(string)
			elseif (imis==84) then
				filename="MIS84"//trim(string)
			elseif (imis==85) then
				filename="MIS85"//trim(string)
			elseif (imis==86) then
				filename="MIS86"//trim(string)
			elseif (imis==87) then
				filename="MIS87"//trim(string)
			elseif (imis==88) then
				filename="MIS88"//trim(string)
			elseif (imis==89) then
				filename="MIS89"//trim(string)
			elseif (imis==90) then
				filename="MIS90"//trim(string)

			elseif (imis==91) then
				filename="MIS91"//trim(string)
			elseif (imis==92) then
				filename="MIS92"//trim(string)
			elseif (imis==93) then
				filename="MIS93"//trim(string)
			elseif (imis==94) then
				filename="MIS94"//trim(string)
			elseif (imis==95) then
				filename="MIS95"//trim(string)
			elseif (imis==96) then
				filename="MIS96"//trim(string)
			elseif (imis==97) then
				filename="MIS97"//trim(string)
			elseif (imis==98) then
				filename="MIS98"//trim(string)
			elseif (imis==99) then
				filename="MIS99"//trim(string)
			elseif (imis==100) then
				filename="MIS100"//trim(string)
			endif
					
			inquire(file=trim(filename), exist=filexist)
			if (.not.filexist) then
				write (*,"(1a)") 'FILE "'//filename//" DOES NOT EXIST, IGNORING THIS PART"
				goto 110
			endif
			nmis_tmp = nmis_tmp+1

!			if (allocated(var)) deallocate(var)
			if (.not.allocated(var)) allocate(var(nmis,nline_max,n))
			write (*,*) "READING FROM "//filename
			open (unit=1,file=trim(filename),status="old",action="read")

!if (iline==1.and.nline==1) read (1,*) zone
			read (1,*) zone

			nline_array(nmis_tmp) = 0

			do iline=1, nline_max
				nline_array(nmis_tmp) = nline_array(nmis_tmp)+1

!				read( 1,*,iostat=io) var(nline_array(nmis_tmp),iline,:)
				read (1,*,iostat=io) var(nmis_tmp, nline_array(nmis_tmp), :)

				if (io>0) then
					write (*,"(1a,1i)") "CHECK "//trim(filename)//". SOMETHING IS WRONG IN LINE ", nline_array(nmis_tmp)
				elseif (io<0) then
					nline_array(nmis_tmp) = nline_array(nmis_tmp) - 1
					write (*,"(1a,1i)") "END OF FILE OCCURED. TOTAL LINES = ", nline_array(nmis_tmp)
					exit
				endif
			enddo



			close (1)
110		continue
		enddo

		if (nmis_tmp<2) goto 10

		if (.not.include_mis) then
			filename = "NMIS"//trim(string)
		else
			filename = "NMIS_"//trim(tmp_str)//trim(string)
		endif
		open (unit=1,file=trim(filename),status="replace",action="write")

		allocate(avr_var(n),var_var(n))

		call get_confin(nmis_tmp, confint)
		write (*,"(1a,1d15.7,1i)") "CONFIDENCE INTERVAL FOR NMIS = ", confint, nmis_tmp

		do iline=1, minval(nline_array(1:nmis_tmp))
			avr_var = 0d0
			var_var = 0d0

			do ivar=1, n
!				avr_var(ivar) = sum(var(:, iline, ivar))

				call calc_avr_var(nmis_tmp, nmis_tmp, var(1:nmis_tmp, iline, ivar), avr_var(ivar), var_var(ivar), .false.)
			enddo
!			avr_var(:) = avr_var(:)/nmis

!			do ivar=1, n
!				do imis=1, nmis
!					var_var(ivar) = var_var(ivar) + (var(imis, iline, ivar)-avr_var(ivar))**2
!				enddo
!			enddo

			if (nmis_tmp>1) then
				var_var(:) = var_var(:) / (nmis_tmp) / (nmis_tmp-1)
				var_var(:) = confint * sqrt(var_var(:))
			else
				var_var(:) = zero
			endif

			if (.not.include_mis) then
				write (1,trim(outform)) avr_var(:), var_var(:)
			else
				write (1,trim(outform)) nmis_tmp, avr_var(:), var_var(:)
			endif
		enddo
		close (1)
10		if (allocated(var)) deallocate(var)
		if (allocated(var_var)) deallocate (avr_var, var_var)
	end subroutine mis_average


	subroutine calc_avr_var(nvar, nvar_total, var, avr, variance, parallel_mode)
		implicit none
		integer, intent(in) :: nvar, nvar_total
		real(prcn), intent(in) :: var(nvar)
		real(prcn), intent(out) :: avr, variance
		logical, intent(in) :: parallel_mode

		integer :: ivar
		real(prcn) :: tmp

		avr = sum(var(:))

		if (parallel_mode) then
		    tmp=avr
		    GLOBAL_DOUBLE_SUM(tmp,avr,1,comm_cart_2d)
		endif
		avr = avr/nvar_total

		variance = zero
		do ivar=1, nvar
			variance = variance + (var(ivar)-avr)**2
		enddo

		if (parallel_mode) then
			tmp=variance
			GLOBAL_DOUBLE_SUM(tmp,variance,1,comm_cart_2d)
		endif
	end subroutine calc_avr_var


	subroutine combine_history
		implicit none
		real(prcn), allocatable, dimension(:,:,:) :: var, var_fin
		real(prcn), allocatable, dimension(:) :: var_tmp
		real(prcn), allocatable, dimension(:,:) :: avr_var, var_var
		integer, allocatable, dimension(:) :: step

		real(prcn) :: time, confint, gtmp1, gtmp2
		integer :: mintime, minstep, maxstep, istep, jstep, imis, ivar, nvar, des_var1, des_var2, ihist, nhist
		integer :: i, j, io, nmis_tmp
		character*50 filename, filename1, filename2, tmp_str
		logical :: filexist, diagnostic=.false.

		integer :: nzones, tmp_nbins, izone, ibin
		integer :: nzones_max = 10000

		real(prcn), allocatable :: time_zone(:), var_zone(:,:,:,:), avr_var_zone(:,:,:), var_var_zone(:,:,:)
	
		nhist = 3
		maxstep = 100000

!if (I_AM_NODE_ZERO) THEN

!		call get_confin(nmis, confint)
!		write (*,*) "CONFIDENCE = ", confint

		allocate(step(nmis))

		do ihist=1, 3
			if (ihist==1) then
				filename="tke.dat"
!				nvar = 24
!				des_var1 = 5
!				des_var2 = 5

!				nvar = 5
!				des_var1 = 1
!				des_var2 = 2

				nvar = 17
				des_var1 = 1
				des_var2 = 1
			elseif (ihist==2) then
				filename="norm_drag.dat"
				if (nphases==1) then
					nvar = 11
				else
					nvar = 12 !4 + nphases + ndim * 2
				endif
				des_var1 = 1
				des_var2 = 1
			elseif (ihist==3) then
				filename="vel_info.dat"
				if (nphases==1) then
					nvar = 12
				else
					nvar = 16 !4+ ndim + ndim * phsc2count + nphases + ndim + 1
				endif
				des_var1 = 1
				des_var2 = 2
			elseif (ihist==4) then
				filename="gofr_vs_t.dat"
!				nvar = 5
				nvar = 3
				des_var1 = 0
				des_var2 = 0
			endif

			nmis_tmp = 0
			do imis=1, nmis
				if (imis==1) then
					filename1="MIS1_"//trim(filename)
				elseif (imis==2) then
					filename1="MIS2_"//trim(filename)
				elseif (imis==3) then
					filename1="MIS3_"//trim(filename)
				elseif (imis==4) then
					filename1="MIS4_"//trim(filename)
				elseif (imis==5) then
					filename1="MIS5_"//trim(filename)
				elseif (imis==6) then
					filename1="MIS6_"//trim(filename)
				elseif (imis==7) then
					filename1="MIS7_"//trim(filename)
				elseif (imis==8) then
					filename1="MIS8_"//trim(filename)
				endif
			
				inquire(file=trim(filename1), exist=filexist)
				if (.not.filexist) then
					write (*,*) 'FILE "'//trim(filename1)//' DOES NOT EXIST, SKIPPING '//trim(filename)
					goto 10
				endif
				nmis_tmp = nmis_tmp+1

				write (*,*) "READING FROM "//filename1
				open (unit=1,file=trim(filename1),status="old",action="read")

				if (ihist==4) then
					nzones = 0
					tmp_nbins = 0

					if (.not.allocated(var_zone)) then
						allocate (time_zone(nzones_max))
						allocate (var_zone(nmis, nzones_max, nvar, nbins))
						allocate (var_tmp(nvar))
						time_zone = zero
						var_zone = zero
						var_tmp = zero
					endif

					do
						read (1,"(1A50)",iostat=io) tmp_str

						if (io>0) then
							write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN ZONE ", nzones+1
						elseif (io==0) then
							filename2 = trim(tmp_str)

!								write (*,*) filename2(2:5), trim(filename2)
!								read (*,*)


							if (filename2(2:5)=="ZONE".or.filename2(2:5)=="zone") then
!write (*,*) "in"
								if (nzones>0.and.tmp_nbins/=nbins) write (*,"(1a,1i5,1a,1i5)") "WARNING! NBINS IN ZONE ",nzones," = ", tmp_nbins
								nzones = nzones + 1
								tmp_nbins = 0
							else
								tmp_nbins = tmp_nbins+1
							endif
						elseif (io<0) then
							write (*,*) "END OF FILE ", trim(filename1)
							exit
						endif
					enddo
					step(imis) = nzones
					write(*,"(1A,2I)")'NUMBER OF ZONES, BINS = ', nzones, nbins
					close(1)
					open (unit=1,file=trim(filename1),status="old",action="read")

					do izone=1, nzones
						read (1,"(1A50)") tmp_str
						do ibin=1, nbins
							read (1,"(3d15.7)",iostat=io) var_tmp(1:nvar)

							if (io>0) then
								write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN ZONE ", izone," LINE ", ibin
							elseif (io==0) then
								var_zone(imis, izone, 1:nvar,ibin) = var_tmp (1:nvar)
							elseif (io<0) then
								write (*,*) "END OF FILE ", trim(filename1)
								exit
							endif
						enddo
						if (io<0) exit
					enddo
				else

					if (.not.allocated(var)) then
						allocate(var(nmis,nvar,maxstep),var_tmp(nvar))
						var = 0
						step = 0
					endif


					step(nmis_tmp) = 0
					do
						if (ihist==1) then
							read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
!							var_tmp(des_var1) = var_tmp(des_var1) / (one-mean_volfrac)
						elseif (ihist==2) then
							if (nphases==1) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
							elseif (nphases==2) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
							endif
							var_tmp(des_var1) = var_tmp(des_var1) / (one-phiavg)
						elseif (ihist==3) then
!							if (nphases==1) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
!write (*,*) var_tmp(8)
								var_tmp(8) = var_tmp(8)**2 * 3 * rhos * phiavg
!write (*,*) var_tmp(8)
!read (*,*)
!							elseif (nphases==2) then
!								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
!								gtmp1 = var_tmp(11)**2 * 3/2 * rhos * phase_array(1)%volfrac
!								gtmp2 = var_tmp(12)**2 * 3/2 * rhos * phase_array(2)%volfrac
!								var_tmp(11) = gtmp1 / (gtmp1+gtmp2)
!								var_tmp(12) = gtmp2 / (gtmp1+gtmp2)
!							endif
!							var_tmp(des_var1) = var_tmp(des_var1) / (one-mean_volfrac)
						endif
						if (io>0) then
							write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN LINE ", step(nmis_tmp)+1
						elseif (io==0) then
							step(nmis_tmp) = step(nmis_tmp) + 1
							var(nmis_tmp,1:nvar,step(nmis_tmp)) = var_tmp (1:nvar)
						elseif (io<0) then
							write (*,*) "END OF FILE ", trim(filename1)
							write (*,*) "NUMBER OF INPUTS = ", step(nmis_tmp)
							exit
						endif
					enddo
					close (1)


					istep=2
					do
						if (var(nmis_tmp,des_var1,istep)<=var(nmis_tmp,des_var1,istep-1).and.&
							var(nmis_tmp,des_var2,istep)<=var(nmis_tmp,des_var2,istep-1)) then
							jstep = istep-1
							do
								if (.not.(var(nmis_tmp,des_var1,istep)<=var(nmis_tmp,des_var1,jstep))) then
									jstep = jstep + 1
									exit
								elseif (jstep==1) then
									exit
								else
									jstep = jstep - 1
								endif
							enddo

							write (*,"(1A,1I,1A,1I,1A,D15.5,1A,D15.7)") "PUTTING ", istep," INTO ", jstep, ", ", var(nmis_tmp,des_var1,istep)," with ", var(nmis_tmp,des_var1,jstep)

							j = jstep
							do i=istep, step(nmis_tmp)
								var(nmis_tmp,1:nvar,j) = var(nmis_tmp,1:nvar,i)
								j = j+1
							enddo
							do i=j, step(nmis_tmp)
								var(nmis_tmp,1:nvar,i) = 0
							enddo
							step(nmis_tmp) = j - 1
							istep = jstep
						endif

						istep = istep +1
						if (istep>step(nmis_tmp)) exit
					enddo
				endif
10				continue
			enddo

			call get_confin(nmis_tmp, confint)
			write (*,"(1a,1d15.7,1i)") "CONFIDENCE INTERVAL FOR NMIS = ", confint, nmis_tmp

			if (ihist<4) then

				mintime = 1
				do imis=2, nmis_tmp
					if (var(imis,des_var1,step(imis))<var(mintime,des_var1,step(mintime))) mintime = imis
				enddo
				minstep = step(mintime)

				write (*,*) "MINIMUM TIME = ", var(mintime,des_var1,step(mintime))
				write (*,*) "NUMBER OF STEPS = ", minstep

				allocate(var_fin(nmis_tmp,nvar,minstep))
				do istep=1, minstep
					var_fin(mintime,1:nvar,istep) = var(mintime,1:nvar,istep)
					time = var(mintime,des_var1,istep)
					do imis=1, nmis_tmp
						if (imis/=mintime) then
							do i=1, step(imis)-1
								if (time>=var(imis,des_var1,i).and.time<=var(imis,des_var1,i+1)) then
									var_fin(imis,1:nvar,istep) =  (var(imis,1:nvar,i+1)-var(imis,1:nvar,i)) / (var(imis,des_var1,i+1)-var(imis,des_var1,i)) &
														&  * (time-var(imis,des_var1,i)) + var(imis,1:nvar,i)
!									var_fin(imis,des_var1,istep) = time
									exit
								endif
							enddo
						endif
					enddo
				enddo
				write (*,*) "INTERPOLATION DONE"

				allocate(avr_var(nvar,minstep),var_var(nvar,minstep))
				avr_var = 0d0
				var_var = 0d0

				do istep=1, minstep
					do ivar=1, nvar
						avr_var(ivar,istep) = sum(var_fin(1:nmis_tmp,ivar,istep))/nmis_tmp
					enddo
				enddo

				write (*,*) "AVERAGE COMPUTED"

				if (nmis_tmp>1) then
					do istep=1, minstep
						do ivar=1, nvar
							do imis=1, nmis_tmp
								var_var(ivar,istep) = var_var(ivar,istep) + (var_fin(imis,ivar,istep)-avr_var(ivar,istep))**2
							enddo
						enddo
						if (nmis_tmp>2) then
							var_var(:,istep) = var_var(:,istep)/nmis_tmp/(nmis_tmp-1)
							var_var(:,istep) = confint*sqrt(var_var(:,istep))
						else
							var_var(:,istep) = var_var(:,istep)/nmis_tmp
							var_var(:,istep) = sqrt(var_var(:,istep))
						endif
					enddo
				else
					var_var(:,istep) = zero
				endif

				write (*,*) "--------------------"

			elseif (ihist==4) then
				minstep = minval(step)
				allocate(avr_var_zone(minstep,nvar,nbins), var_var_zone(minstep,nvar,nbins))
				avr_var_zone = 0d0
				var_var_zone = 0d0

				do izone=1, minstep
					do ibin=1, nbins
						do ivar=1, nvar
							avr_var_zone(izone, ivar,ibin) = sum(var_zone(1:nmis, izone, ivar, ibin)) / nmis
						enddo
					enddo


					if (nmis>1) then
						do ibin=1, nbins
							do ivar=1, nvar
								do imis=1, nmis
									var_var_zone(izone, ivar, ibin) = var_var_zone(izone, ivar, ibin) + &
									&(var_zone(imis, izone, ivar, ibin)-avr_var_zone(izone, ivar, ibin))**2
								enddo
							enddo
							var_var_zone(izone, : , ibin) = var_var_zone(izone, : , ibin) / nmis/(nmis-1)
							var_var_zone(izone, : , ibin) = confint * sqrt(var_var_zone(izone, : , ibin))
						enddo
					else
						var_var(:,istep) = zero
					endif

				enddo

				write (*,*) "AVERAGE COMPUTED"
				write (*,*) "95% CONFIDENCE INTERVAL COMPUTED"
				write (*,*) "--------------------"
			endif

			filename2 = "NMIS_"//trim(filename)
			open (unit=1,file=trim(filename2),status="replace",action="write")

			if (ihist<4) then
				do istep=1, minstep
					if (ihist==1) then
						write (1,"(2(24D15.7))") avr_var(:,istep), var_var(:,istep)
					elseif (ihist==2) then
						write (1,"(24(2x,e20.12))") avr_var(:,istep), var_var(:,istep)
					elseif (ihist==3) then
						write (1,"(32(2x,g17.8))") avr_var(:,istep), var_var(:,istep)
					endif
				enddo
			else
				do izone=1, minstep
					write (1,*) "ZONE"
					do ibin=1, nbins
						if (ihist==4) then
							write (1,"(6D15.7)") avr_var_zone(izone,:,ibin), var_var_zone(izone,:,ibin)
						endif
					enddo
				enddo
			endif

			close (1)
20			if (allocated(var)) deallocate(var, var_tmp, var_fin, avr_var, var_var)
			if (allocated(var_zone)) deallocate(time_zone, var_zone, var_tmp, avr_var_zone, var_var_zone)
		enddo
		deallocate(step)
!endif
	end subroutine combine_history


	SUBROUTINE flow_snapshot
		use nlmainarrays, only : ubcp, pbcp
		USE dem_mod, only : is_mobile, des_pos_new, des_radius
		IMPLICIT NONE
		Integer  :: sunit,i,j,k,l,m,isp, mark, idim
		INTEGER, SAVE :: zone_count = 0
		LOGICAL, SAVE :: first_pass=.TRUE.
		REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
			& mean_vel(ndim), mean_force(ndim), position(ndim)
		CHARaCTER*100 :: FILENAME1, filename2, filename3
		integer, save :: sphrunit1, sphrunit2, sphrunit3

		real(prcn), allocatable :: out_arr(:,:,:,:), trans_buf(:)
		integer :: node_num, iproc
		integer :: ii, jj, kk, j1, j2, j3
		real(prcn) :: tmp, tmp1, tmp2, tmp3, mean_energy
		integer :: iphs, part_start, part_end
		logical :: filexist

		real(prcn), allocatable :: acc_fluc(:,:), vel_fluc(:,:), acc_var(:,:), vel_var(:,:), acc_fluc_meanf(:,:), acc_var_meanf(:,:)

		j1=1
		j2=global_n(2)/2
		j3=global_n(2)

		if (I_AM_NODE_ZERO) write (*,*) "IN FLOW_SNAPSHOT"

#if PARALLEL
		if (I_AM_NODE_ZERO) then
			allocate(out_arr(global_n(1), 3, global_n(3), ndim+3))

			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				! velocity fluctuations for node zero
				if (local_i_start(2)<=jj .and. jj<=local_i_end(2)) then
					do idim=1, ndim+1
						do i=1, local_ni(1)
							do k=1, local_ni(3)

								ii = local_i_start(1) + i - 1
								kk = local_i_start(3) + k - 1

								if (idim<=ndim) then
									out_arr(ii,jj,kk,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
									!if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
								else
									out_arr(ii,jj,kk,ndim+1) = pbcp(i,j,k) !-ufmean(idim)
								endif
							enddo
						enddo
					enddo
				endif
			enddo

			do iproc=1, nproc-1
				if (iystarts(iproc)<=jj .and. jj<=iyends(iproc)) then
					node_num = izlocals(iproc) * ixlocals(iproc) * (ndim+1)
					allocate(trans_buf(node_num))

					do jj=1, 3
						if (jj==1) then
							j=j1
						elseif (jj==2) then
							j = j2
						elseif (jj==3) then
							j = j3
						endif

						! collecting velocity fluctuations from other processes
						call mpi_recv(trans_buf(1), node_num, mpi_double_precision, iproc, iproc, comm_cart_2d, status, err_code)

						l=0
						do idim=1, ndim+1
							do k=izstarts(iproc), izends(iproc)
								do i=ixstarts(iproc),ixends(iproc)
									l=l+1
									if (idim<=ndim) then
										out_arr(i,jj,k,idim) = trans_buf(l)
									else
										out_arr(i,jj,k,ndim+1) = trans_buf(l)
									endif
								enddo
							enddo
						enddo
					enddo
					deallocate(trans_buf)
				endif
			enddo
		else
			! recieving velocity fluctuations from node zero

			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				if (local_i_start(2)<=jj .and. jj<=local_i_end(2)) then

					node_num=local_ni(3)*local_ni(1)*(ndim+1)
					allocate(trans_buf(node_num))

					l=0
					do idim=1,ndim+1
						do k=1, local_ni(3)
							do i=1, local_ni(1)
								l=l+1
								if (idim<=ndim) then
									trans_buf(l) = ubcp(i,j,k,idim) !-ufmean(idim)
									!if (.not.fluid_atijk(i,j,k)) trans_buf(l) = 0d0
								else
									trans_buf(l) = pbcp(i,j,k) !-ufmean(idim)
								endif
							enddo
						enddo
					enddo

					call mpi_send(trans_buf(1), node_num, mpi_double_precision, node_zero, myid, comm_cart_2d, err_code)

					deallocate(trans_buf)
				endif
			enddo

		endif
#else
		allocate(out_arr(global_n(1),3,global_n(3),ndim+3))
		do idim=1, ndim+1
			do k=1, global_n(3)
				do jj=1, 3
					if (jj==1) then
						j=j1
					elseif (jj==2) then
						j = j2
					elseif (jj==3) then
						j = j3
					endif

					do i=1, global_n(1)
						if (idim<=ndim) then
							out_arr(i,jj,k,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
	!						if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
						else
							out_arr(i,jj,k,idim) = pbcp(i,j,k)
						endif
					enddo
				enddo
			enddo
		enddo
#endif
		mean_energy = 0.5*umeanslip**2
		write (*,*) "mean_slip, mean_energy = ", umeanslip, mean_energy
		if (I_AM_NODE_ZERO) then
			write (*,*) "GENERATING THE SNAPSHOT OF THE FIELD"

			do k=1, global_n(3)
				do j=1, 3
					do i=1, global_n(1)
!						out_arr(i,j,k,1:ndim) = out_arr(i,j,k,1:ndim) / umeanslip
!						out_arr(i,j,k,ndim+1) = out_arr(i,j,k,ndim+1) / mean_energy
						out_arr(i,j,k,ndim+3) = abs(dot_product(out_arr(i,j,k,1:ndim),out_arr(i,j,k,1:ndim)))

						tmp = 0d0
						do idim=1, ndim
							tmp = tmp + (out_arr(i,j,k,idim)-ufmean(idim)) * (out_arr(i,j,k,idim)-ufmean(idim))
						enddo
						out_arr(i,j,k,ndim+2) = tmp / 2

						out_arr(i,j,k,ndim+3) = out_arr(i,j,k,ndim+1) + mp(1) * (i-global_n(1)/2) /global_n(1)* doml(1)
					enddo
				enddo
			enddo

			FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion_pas'
			filename2 = TRIM(RUN_NAME)//'_sphr_motion_act1'
			filename3 = TRIM(RUN_NAME)//'_sphr_motion_act2'

			sunit     = 30
			filename1 = TRIM(RUN_NAME)//'_MOVIE'
			call instant_file_opener(filename1, sunit, .true.)

			write(sunit,*)'ZONE T = "', t/t_conv, '",'
			write(sunit,"(3(1a,1i))")'DATAPACKING=POINT, I =', global_n(1),  ', J=', 3, ', K=', global_n(3)

			do k=1, global_n(3)
				do jj=1, 3
					if (jj==1) then
						j=j1
					elseif (jj==2) then
						j = j2
					elseif (jj==3) then
						j = j3
					endif

					do i=1, global_n(1)
						write(sunit,"(3i,10d15.7)") i, j, k, out_arr(i,jj,k,1:ndim) / umeanslip , out_arr(i,jj,k,4), out_arr(i,jj,k,6) / (vis*umeanslip/dia_phys)
					enddo
				enddo
			enddo
			close(sunit)

			sphrunit1 = 31
			filename2 = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'
			call instant_file_opener(filename2, sphrunit1, .true.)

			WRITE(sphrunit1,*)'ZONE T= "', t/t_conv, ' " '

			do iphs=1, nphases
				part_start = phase_array(iphs)%pstart
				part_end   = phase_array(iphs)%pend

				do m=part_start, part_end
					tmp1 = sqrt(dot_product(velbdy(m,:),velbdy(m,:))) / umeanslip
					!tmp2 = dot_product(velbdy(m,:)-phase_array(iphs)%mean_spec_vel(:), velbdy(m,:)-phase_array(iphs)%mean_spec_vel(:)) / 2 / mean_energy
					!if (dot_product(acc_var(iphs,:),vel_var(iphs,:))>post_small) then
					!	tmp3 = dot_product(acc_fluc(m,:),vel_fluc(m,:)) / dot_product(acc_var(iphs,:),vel_var(iphs,:))
					!else
					!	tmp3 = dot_product(acc_fluc(m,:),vel_fluc(m,:))
					!endif

					write (sphrunit1,"(6d15.7)")  xc(m,:), radbdy(m), tmp1 !, tmp2
!					write (sphrunit1,"(7d15.7)")  xc(m,:), radbdy(m), tmp1, tmp2, tmp3
				enddo
			enddo
			close(sphrunit1)

			deallocate(out_arr)
		endif
	END SUBROUTINE flow_snapshot








	subroutine uiui_correlation
		use postproc_funcs
		!use init_turb, only : gener_filename
		implicit none
		integer :: i, j, k, l, i1, i2, j1, j2, k1, k2, ii1, ii2, imis, idim, nvar, nvar1, nvar2, nline, ivar, ibin
		integer(8) :: ijk1, ijk2, ijk_start
		real(prcn) :: dr, rmax, var_lint, avr_lint, sd_lint, err_bar, confint, rij_mag
		real(prcn), allocatable :: uiui(:,:), uiui_par(:), uiui_per(:)
		real(prcn), allocatable :: uiui_ss(:,:), uiui_ss_par(:), uiui_ss_per(:)
		real(prcn), allocatable :: uiui_fs(:,:), uiui_fs_par(:), uiui_fs_per(:)
		real(prcn), allocatable :: uiui_sf(:,:), uiui_sf_par(:), uiui_sf_per(:)

		integer(8), allocatable :: num_bins(:), num_bins_fs(:), num_bins_sf(:), num_bins_ss(:)
		real(prcn) :: ui(ndim), uj(ndim), u1i(ndim), u1j(ndim), u2i(ndim), u2j(ndim), u1i_mag, u1j_mag
		real(prcn) :: cpu1, cpu0

		integer :: iproc, id_s, id_r, node_num_s, node_num_r, s_size
		real(prcn), allocatable, dimension(:,:,:,:) :: urecv

		logical, allocatable :: fluid_atijk2(:,:,:)


#if PARALLEL
		integer, allocatable :: samplearray(:)
		integer :: ierr
		integer :: mysample
#endif

		character*50 filename1
		character*5 turn
		character*1 rank_string
		integer :: proc_start, count_out, count=0
		integer :: strlen
		logical :: filexist, finish

		rmax  = doml(2) / 2

		allocate(uiui(nbins,ndim), uiui_par(nbins), uiui_per(nbins))
		allocate(num_bins(nbins))

		uiui     = zero
		uiui_par = zero
		uiui_per = zero
		num_bins = 0
		


		if (imove==1.or.move_particles) then
			allocate(uiui_ss(nbins,ndim), uiui_ss_par(nbins), uiui_ss_per(nbins))
			allocate(uiui_fs(nbins,ndim), uiui_fs_par(nbins), uiui_fs_per(nbins))
			allocate(uiui_sf(nbins,ndim), uiui_sf_par(nbins), uiui_sf_per(nbins))

			uiui_ss     = zero
			uiui_ss_par = zero
			uiui_ss_per = zero
			num_bins_ss = 0

			uiui_fs     = zero
			uiui_fs_par = zero
			uiui_fs_per = zero
			num_bins_fs = 0

			uiui_sf     = zero
			uiui_sf_par = zero
			uiui_sf_per = zero
			num_bins_sf = 0
		endif

		if (allocated(rad_bin)) deallocate(rad_bin)
		allocate(rad_bin(nbins))
		rad_bin = zero

		dr = rmax / nbins

		id_r = myid
		id_s = myid

		do ibin=1, nbins
			rad_bin(ibin) = (ibin-0.5) * dr
		enddo

#if PARALLEL

		allocate(samplearray(0:nproc-1))

		mysample = local_ni(3)*local_ni(2)*local_ni(1)

		call mpi_allgather(mysample, 1, mpi_int, samplearray, 1, mpi_int, comm_cart_2d, ierr)
#endif

		if (.not.post_no_flow_mem_alloc) then
			call screen_separator(30,'I')
			if (I_AM_NODE_ZERO) write (*,*) "GENERATING THE PARALLEL AND PERPENDICULAR VELOCITY CORRELATIONS"

			call initialize_gridvertindex(local_ni(1), local_ni(2), local_ni(3))

			call GENER_FILENAME(filename1, trim(run_name)//"_corr_res.rst")
			inquire(file=filename1, exist=filexist)

			if (filexist) then
				call restart_in
			else
				proc_start = 0
				ijk_start  = 0
			endif

			do iproc = proc_start, nproc/2
				count_out = 0
				first_pass = .true.
				do ijk1 = ijk_start+1, nvert, skip_num
					count_out = count_out+1
					call ind1t3(ijk1,i1,j1,k1)

					if (.not.((imove==1.or.move_particles) .or. fluid_atijk(i1,j1,k1))) goto 20

					if (iproc==0) then
						do ijk2 = ijk1, nvert
							call ind1t3(ijk2,i2,j2,k2)						
							!if (fluid_atijk(i2,j2,k2)) call calc_correlation(ubcp(1:nx,:,:,:), ubcp(1:nx,:,:,:), i1, j1, k1, i2, j2, k2)
							if ((imove==1 .or. move_particles) .or. fluid_atijk(i2,j2,k2)) then
								call calc_correlation(ubcp(i1,j1,k1,1:ndim), ubcp(i2,j2,k2,1:ndim), fluid_atijk(i1,j1,k1), fluid_atijk(i2,j2,k2), i1, j1, k1, i2, j2, k2)
							endif
						enddo
#if PARALLEL
					else
						if (first_pass) then
							id_s = myid+iproc
							id_r = myid-iproc

							if (id_s>nproc-1) id_s = id_s-nproc
							if (id_r<0)       id_r = id_r+nproc

							node_num_s = mysample
							node_num_r = samplearray(id_r)

							if (allocated(urecv)) deallocate(urecv, fluid_atijk2)
							allocate(urecv(ixstarts(id_r):ixends(id_r), iystarts(id_r):iyends(id_r), izstarts(id_r):izends(id_r), ndim))
							allocate(fluid_atijk2(ixstarts(id_r):ixends(id_r), iystarts(id_r):iyends(id_r), izstarts(id_r):izends(id_r)))

							CALL MPI_SENDRECV(ubcp(1:local_ni(1), 1:local_ni(2), 1:local_ni(3), 1:ndim), node_num_s*ndim, MPI_DOUBLE_PRECISION, id_s, myid, urecv, node_num_r*ndim, MPI_DOUBLE_PRECISION, id_r, id_r, comm_cart_2d, status, err_code)

							CALL MPI_SENDRECV(fluid_atijk(1:local_ni(1), 1:local_ni(2), 1:local_ni(3)), node_num_s, MPI_LOGICAL, id_s, id_s, fluid_atijk2, node_num_r, MPI_LOGICAL, id_r, myid, comm_cart_2d, status, err_code)

							first_pass = .false.
						endif

						if (iproc==nproc/2 .and. myid+1>nproc/2) goto 10

						do k2=izstarts(id_r), izends(id_r)
							do j2=iystarts(id_r), iyends(id_r)
								do i2=ixstarts(id_r), ixends(id_r)
									!if (fluid_atijk2(i2,j2,k2)) call calc_correlation(ubcp(1:nx,:,:,:), urecv, starts(myid)+i1-1, j1, k1, i2, j2, k2)
									if ((imove==1.or.move_particles) .or. fluid_atijk(i2,j2,k2)) then
										call calc_correlation(ubcp(i1,j1,k1,1:ndim), urecv(i2,j2,k2,1:ndim), fluid_atijk(i1,j1,k1), fluid_atijk2(i2,j2,k2), ixstarts(myid)+i1-1, iystarts(myid)+j1-1, izstarts(myid)+k1-1,  i2, j2, k2)
									endif
								enddo
							enddo
						enddo

10						continue
#endif
					endif

20					continue

!					if (mod(count_out,10) == 0) then
!						if (I_AM_NODE_ZERO) write (*,*) "WRITING OUTPUT AND RESTART FILE"
!
!						call output
!						call restart_out
!					endif

					if (I_AM_NODE_ZERO.and.mod(count_out,20)==1) then
#if PARALLEL
						write (*,"(4(1a,1i))") " NODE ", ijk1, " OUT OF ", local_ni(1)*local_ni(2)*local_ni(3), " BETWEEN PROCs ", myid, " AND ", id_r
#else
						write (*,"(2(1a,1i))") " NODE ", ijk1, " OUT OF ", local_ni(1)*local_ni(2)*local_ni(3)
#endif

						call output
						call restart_out

						if (I_AM_NODE_ZERO) write (*,*) "WRITING OUTPUT AND RESTART FILE"
					endif

					if (I_AM_NODE_ZERO) then
						call cpu_time(cpu1)
						cpu1 = cpu1-cpu0
						cpu1 = cpu1/3600
					endif
					BROADCAST_DOUBLE(cpu1, 1, NODE_ZERO, comm_cart_2d)

					killjob = .false.
					if (cpu1>=WTIME_MAXHRS) then
						killjob = .true.
						exit
					endif
				enddo
				if (killjob) exit
				ijk_start = 0
			enddo

			if (killjob) then
				if (I_AM_NODE_ZERO) then
					write (*,*) "RUNNING EXCEEDED THE WTIME_MAXHRS"
					write (*,*) "CORRELATION COMPATION NOT FINISHED YET"
					write (*,*) "RESTART THE JOB AGAIN"
				endif
			endif

30			call output
			call restart_out

			deallocate (uiui, uiui_par, uiui_per)
			deallocate (rad_bin, num_bins)

			if (imove==1.or.move_particles) then
				deallocate (uiui_ss, uiui_ss_par, uiui_ss_per, num_bins_ss)
				deallocate (uiui_fs, uiui_fs_par, uiui_fs_per, num_bins_fs)
				deallocate (uiui_sf, uiui_sf_par, uiui_sf_per, num_bins_sf)
			endif

			if (I_AM_NODE_ZERO) write (*,"(1A)") "EXITING UIUI_CORRELATION"
		else
			if(I_AM_NODE_ZERO)then
				if (imove==1.or.move_particles) then
					nvar  = 13
					nvar1 = 1
					nvar2 = 12
				else
					nvar  = 5
					nvar1 = 1
					nvar2 = 4
				endif

				filename1 = "_uiui_correlation.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line, include_nmis)
			endif 
		endif

		if (I_AM_NODE_ZERO) call screen_separator(30,'I')

	contains

		subroutine calc_correlation(u1, u2, fluid1, fluid2, i1, j1, k1, i2, j2, k2)
			implicit none
			integer, intent(in) :: i1, i2, j1, j2, k1, k2
			real(prcn), intent(in) :: u1(ndim), u2(ndim)
			logical, intent(in) :: fluid1, fluid2
			real(prcn) :: r, rij(ndim), vec1(ndim), vec2(ndim), vec1_par(ndim), vec2_par(ndim), vec1_per(ndim), vec2_per(ndim), unit(ndim)
			integer :: body1, body2

			if (.not.fluid1) call find_body(i1,j1,k1,body1, rmax)
			if (.not.fluid2) call find_body(i2,j2,k2,body2, rmax)

			rij(1) = (i2-i1) * dx
			rij(2) = (j2-j1) * dy
			rij(3) = (k2-k1) * dz


			do idim=1, ndim
				if (rij(idim) > rmax) then
					rij(idim) = rij(idim)-2*rmax
				elseif (rij(idim) < -rmax) then
					rij(idim) = rij(idim)+2*rmax
				endif
			enddo

			r = sqrt(rij(1)**2 + rij(2)**2 + rij(3)**2) 
			if (r<=rmax) then
				ibin = int(r/rmax * nbins) +1
				if (ibin>nbins) ibin = nbins
				!if (i1==i2.and.j1==j2.and.k1==k2) ibin = 0

				if (fluid1) then
					vec1(:) = u1(:)-ufmean(:)
				else
					vec1(:) = velbdy(body1,:) - usmean(:)
				endif

				if (fluid2) then
					vec2(:) = u2(:)-ufmean(:)
				else
					vec2(:) = velbdy(body2,:) - usmean(:)
				endif

				!if (ibin/=0) then
				if (i1/=i2.or.j1/=j2.or.k1/=k2) then
					unit(:) = rij(:) / r
					vec1_par(:) = dot_product(vec1(:),unit(:)) * unit(:)
					vec1_per(:) = vec1(:) - vec1_par(:) 
					vec2_par(:) = dot_product(vec2(:),unit(:)) * unit(:)
					vec2_per(:) = vec2(:) - vec2_par(:) 
				else
					vec1_par = zero
					vec1_per = zero
					vec2_par = zero
					vec2_per = zero

					vec1_par(1) = vec1(1)
					vec1_per(2) = vec1(2)
					vec2_par(1) = vec2(1)
					vec2_per(2) = vec2(2)
				endif

				if (fluid1.and.fluid2) then
					num_bins(ibin)  = num_bins(ibin) + 1

					uiui(ibin,:)   = uiui(ibin,:)   + vec1(:) * vec2(:)
					uiui_par(ibin) = uiui_par(ibin) + dot_product(vec1_par, vec2_par)
					uiui_per(ibin) = uiui_per(ibin) + dot_product(vec1_per, vec2_per)
				elseif ((.not.fluid1) .and. (.not.fluid2)) then
					num_bins_ss(ibin)  = num_bins_ss(ibin) + 1

					uiui_ss(ibin,:)   = uiui_ss(ibin,:)   + vec1(:) * vec2(:)
					uiui_ss_par(ibin) = uiui_ss_par(ibin) + dot_product(vec1_par, vec2_par)
					uiui_ss_per(ibin) = uiui_ss_per(ibin) + dot_product(vec1_per, vec2_per)
				elseif (fluid1 .and. (.not.fluid2)) then
					num_bins_fs(ibin)  = num_bins_fs(ibin) + 1

					uiui_fs(ibin,:)   = uiui_fs(ibin,:)   + vec1(:) * vec2(:)
					uiui_fs_par(ibin) = uiui_fs_par(ibin) + dot_product(vec1_par, vec2_par)
					uiui_fs_per(ibin) = uiui_fs_per(ibin) + dot_product(vec1_per, vec2_per)
				elseif ((.not.fluid1) .and. fluid2) then
					num_bins_sf(ibin)  = num_bins_sf(ibin) + 1

					uiui_sf(ibin,:)   = uiui_sf(ibin,:)   + vec1(:) * vec2(:)
					uiui_sf_par(ibin) = uiui_sf_par(ibin) + dot_product(vec1_par, vec2_par)
					uiui_sf_per(ibin) = uiui_sf_per(ibin) + dot_product(vec1_per, vec2_per)
				endif
			endif
		end subroutine calc_correlation

		subroutine output
			implicit none
			character*50 filename
			real(prcn), dimension(nbins,ndim) :: u1out, u1out_ss, u1out_fs, u1out_sf
			real(prcn), dimension(nbins) :: u2out, u3out, u2out_ss, u3out_ss, u2out_fs, u3out_fs, u2out_sf, u3out_sf
			integer(8), dimension(nbins) :: num_binsout, num_binsout_fs, num_binsout_sf, num_binsout_ss
			integer :: i, j

			GLOBAL_DOUBLE_SUM(uiui, u1out, (nbins)*ndim, comm_cart_2d)
			GLOBAL_DOUBLE_SUM(uiui_par, u2out, nbins, comm_cart_2d)
			GLOBAL_DOUBLE_SUM(uiui_per, u3out, nbins, comm_cart_2d)
			GLOBAL_LONGINT_SUM(num_bins, num_binsout, nbins, comm_cart_2d)

			if (imove==1.or.move_particles) then
				GLOBAL_DOUBLE_SUM(uiui_ss, u1out_ss, (nbins)*ndim, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(uiui_ss_par, u2out_ss, nbins, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(uiui_ss_per, u3out_ss, nbins, comm_cart_2d)
				GLOBAL_LONGINT_SUM(num_bins_ss, num_binsout_ss, nbins, comm_cart_2d)

				GLOBAL_DOUBLE_SUM(uiui_fs, u1out_fs, (nbins)*ndim, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(uiui_fs_par, u2out_fs, nbins, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(uiui_fs_per, u3out_fs, nbins, comm_cart_2d)
				GLOBAL_LONGINT_SUM(num_bins_fs, num_binsout_fs, nbins, comm_cart_2d)

				GLOBAL_DOUBLE_SUM(uiui_sf, u1out_sf, (nbins)*ndim, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(uiui_sf_par, u2out_sf, nbins, comm_cart_2d)
				GLOBAL_DOUBLE_SUM(uiui_sf_per, u3out_sf, nbins, comm_cart_2d)
				GLOBAL_LONGINT_SUM(num_bins_sf, num_binsout_sf, nbins, comm_cart_2d)
			endif

			if (I_AM_NODE_ZERO) then
				do ibin=1, nbins
					if (num_binsout(ibin)>0) then
						u1out(ibin,:) = u1out(ibin,:) / num_binsout(ibin)
						u2out(ibin) = u2out(ibin) / num_binsout(ibin)
						u3out(ibin) = u3out(ibin) / num_binsout(ibin)
					endif

					if (imove==1.or.move_particles) then
						if (num_binsout_ss(ibin)>0) then
							u1out_ss(ibin,:) = u1out_ss(ibin,:) / num_binsout_ss(ibin)
							u2out_ss(ibin) = u2out_ss(ibin) / num_binsout_ss(ibin)
							u3out_ss(ibin) = u3out_ss(ibin) / num_binsout_ss(ibin)
						endif

						if (num_binsout_fs(ibin)>0) then
							u1out_fs(ibin,:) = u1out_fs(ibin,:) / num_binsout_fs(ibin)
							u2out_fs(ibin) = u2out_fs(ibin) / num_binsout_fs(ibin)
							u3out_fs(ibin) = u3out_fs(ibin) / num_binsout_fs(ibin)
						endif

						if (num_binsout_sf(ibin)>0) then
							u1out_sf(ibin,:) = u1out_sf(ibin,:) / num_binsout_sf(ibin)
							u2out_sf(ibin) = u2out_sf(ibin) / num_binsout_sf(ibin)
							u3out_sf(ibin) = u3out_sf(ibin) / num_binsout_sf(ibin)
						endif
					endif
				enddo

				do ibin=nbins ,1, -1
					do idim=1, ndim
						if (u1out(1,idim)>small_number) u1out(ibin,idim) = u1out(ibin,idim) / u1out(1,idim)
					enddo
					if (u2out(1)>small_number) u2out(ibin) = u2out(ibin) / u2out(1)
					if (u3out(1)>small_number) u3out(ibin) = u3out(ibin) / u3out(1)

					if (imove==1.or.move_particles) then
						do idim=1, ndim
							if (u1out_ss(1,idim)>small_number) u1out_ss(ibin,idim) = u1out_ss(ibin,idim) / u1out_ss(1,idim)
						enddo
						if (u2out_ss(1)>small_number) u2out_ss(ibin) = u2out_ss(ibin) / u2out_ss(1)
						if (u3out_ss(1)>small_number) u3out_ss(ibin) = u3out_ss(ibin) / u3out_ss(1)

						do idim=1, ndim
							if (u1out_fs(1,idim)>small_number) u1out_fs(ibin,idim) = u1out_fs(ibin,idim) / u1out_fs(1,idim)
						enddo
						if (u2out_fs(1)>small_number) u2out_fs(ibin) = u2out_fs(ibin) / u2out_fs(1)
						if (u3out_fs(1)>small_number) u3out_fs(ibin) = u3out_fs(ibin) / u3out_fs(1)

						do idim=1, ndim
							if (u1out_sf(1,idim)>small_number) u1out_sf(ibin,idim) = u1out_sf(ibin,idim) / u1out_sf(1,idim)
						enddo
						if (u2out_sf(1)>small_number) u2out_sf(ibin) = u2out_sf(ibin) / u2out_sf(1)
						if (u3out_sf(1)>small_number) u3out_sf(ibin) = u3out_sf(ibin) / u3out_sf(1)
					endif
				enddo

				filename = trim(run_name)//"_uiui_correlation.dat"

!				inquire(file=trim(filename), exist=filexist)
!				if (filexist) then
!					open(unit=1, file=trim(filename), status="old", position="append")
!				else
					open(unit=1, file=trim(filename), status="replace")
!				endif
				write (1,*) "zone"

				do ibin=1, nbins
					if (imove==1.or.move_particles) then
						write (1,'(1f10.4,20d15.7)') rad_bin(ibin), u1out(ibin,:), u2out(ibin), u3out(ibin), &
														&                   u1out_ss(ibin,:), u2out_ss(ibin), u3out_ss(ibin), &
														&                   u1out_fs(ibin,:), u2out_fs(ibin), u3out_fs(ibin), &
														&                   u1out_sf(ibin,:), u2out_sf(ibin), u3out_sf(ibin)
					else
						write (1,'(1f10.4,5d15.7)') rad_bin(ibin), u1out(ibin,:), u2out(ibin), u3out(ibin)
					endif
				enddo
				close (1)
			endif
		end subroutine output

		subroutine restart_in
			implicit none

			character*50 filename1, filename2

			call GENER_FILENAME(filename1,trim(run_name)//"_corr_res.rst")

			open (unit=9998,file=filename1,status="old",action="read")
			read (9998,*) count
			close (9998)

			if (count==0) then
				turn = "_0"
			else
				turn = "_1"
			endif
#if PARALLEL
			if (I_AM_NODE_ZERO) then
				do i=nproc-1,0,-1
					filename2 = trim(run_name)//"_corr_res"//trim(turn)
					call to_string(i, rank_string, s_size)
!						do j=1, 3-s_size
!							rank_string="0"//trim(rank_string)
!						enddo
					filename2 = trim(filename2)//"_"//trim(rank_string)//".rst"
					if (i/=node_zero) SEND_STRING(filename2, strlen, i, 0, 1, comm_cart_2d)
				enddo
			else
				RECV_STRING(filename2, strlen, node_zero, 0, 1, comm_cart_2d, status)
			endif
#else
			filename2 = trim(run_name)//"_corr_res"//trim(turn)//".rst"
#endif	

			if (I_AM_NODE_ZERO) write (*,*) "RESTARTING GENERATION OF AUTOCORRELATION FUNCTION"

			open (unit=9996,file=filename2,status="old",action="read",form="unformatted")

			if (imove==1.or.move_particles) then
				read (9996) proc_start, ijk_start, uiui, uiui_par, uiui_per, num_bins, uiui_ss, uiui_ss_par, uiui_ss_per, num_bins_ss, uiui_fs, uiui_fs_par, uiui_fs_per, num_bins_fs, uiui_sf, uiui_sf_par, uiui_sf_per, num_bins_sf
 
			else
				read (9996) proc_start, ijk_start, uiui, uiui_par, uiui_per, num_bins
			endif
			close (9996)
		end subroutine restart_in

		subroutine restart_out
			implicit none
			character*50 filename
			character*50 filename1

			count = count+1

			if (mod(count,2)==0) then
				turn = "_0"
			else
				turn = "_1"
			endif

#if PARALLEL
			if (I_AM_NODE_ZERO) then
				do i=nproc-1,0,-1
					filename = trim(run_name)//"_corr_res"//trim(turn)
					call to_string(i, rank_string, s_size)
!						do j=1, 3-s_size
!							rank_string="0"//trim(rank_string)
!						enddo
					filename = trim(filename)//"_"//trim(rank_string)//".rst"
					if (i/=node_zero) SEND_STRING(filename, strlen, i, 0, 1, comm_cart_2d)
				enddo
			else
				RECV_STRING(filename, strlen, node_zero, 0, 1, comm_cart_2d, status)
			endif
#else
			filename = trim(run_name)//"_corr_res"//trim(turn)//".rst"
#endif
			open (unit=9998,file=trim(filename),status="replace",action="write",form="unformatted")
			if (imove==1.or.move_particles) then
				write (9998) iproc, ijk1, uiui, uiui_par, uiui_per, num_bins, uiui_ss, uiui_ss_par, uiui_ss_per, num_bins_ss, uiui_fs, uiui_fs_par, uiui_fs_per, num_bins_fs, uiui_sf, uiui_sf_par, uiui_sf_per, num_bins_sf
			else
				write (9998) iproc, ijk1, uiui, uiui_par, uiui_per, num_bins
			endif
			close (9996)
			close (9998)

			if (I_AM_NODE_ZERO) then
				filename1 = trim(run_name)//"_corr_res.rst"
				open (unit=1,file=trim(filename1),status="replace",action="write")
				write (1,*) mod(count,2)
				close (1)
			endif
		end subroutine restart_out
	end subroutine uiui_correlation


	subroutine find_body(i,j,k,m, rmax)
		implicit none
		integer, intent(in)  :: i, j, k
		real(prcn), intent(in) :: rmax
		integer, intent(out) :: m
		integer :: ibody, idim
		real(prcn) :: dist(ndim), r

		m = 0
		do ibody=1, nbody

			dist(1) = abs(xc(ibody,1)-i) * dx
			dist(2) = abs(xc(ibody,2)-j) * dy
			dist(3) = abs(xc(ibody,3)-k) * dz

			do idim=1, ndim
				if (dist(idim)>rmax) dist(idim) = 2*rmax - dist(idim)
			enddo

			r = sqrt( dot_product(dist(:),dist(:)) )

			if (r<=radbdy(ibody)/dbydx) then
				m = ibody
				exit
			endif
		enddo

		if (m==0) then
			write (*,"(2i,1a)") myid, "AN ERROR IN FINDING THE A SOLID POINT IN A PARTICLE"
		endif
	end subroutine find_body




	subroutine pressure_velocity_vs_theta
		use dependent_functions , only : interpolate_pdata, interpolate_udata
		use bcsetarrays, only :  omega => fr
		use boundary_condition, only : compute_omega
		implicit none

		integer :: l, m, n, iphs, vcellb(ndim), pcellb(ndim), index_out(ndim), ib, ie, jb, je, kb, ke, onew, ii, jj, kk, i, j, k
		real(prcn) :: da, rad, xl(ndim), position_out(ndim), rad2, ul(ndim), ppll(ndim), pl, nll(ndim), onll(ndim), dfll(ndim)
		logical :: i_have_the_point

		real(prcn), allocatable :: theta_bin(:)
		real(prcn), allocatable :: ptheta_local(:), utheta_local(:,:), prestheta_local(:,:), visctheta_local(:,:)
		integer(8), allocatable :: theta_count_local(:)

		real(prcn), allocatable :: ptheta(:), utheta(:,:), prestheta(:,:), visctheta(:,:)
		integer(8), allocatable :: theta_count(:)

		real(prcn) :: dtheta, thetamin, thetamax, xtemp, ytemp, ztemp, phi, theta
		real(prcn) :: df(nbnd,ndim), vort(ndim)
		integer :: itheta, unitnum
		character*100 filename

		if (I_AM_NODE_ZERO) then
			call screen_separator(80,'T')
			write (*,*) "IN PU_vs_THERA..."
		endif


		thetamin = 0
		thetamax = pi
		dtheta = (thetamax-thetamin)/nbins

		allocate(theta_bin(nbins))
		allocate(theta_count_local(nbins))
		allocate(theta_count(nbins))

		allocate(ptheta_local(nbins), utheta_local(nbins,ndim))
		allocate(ptheta(nbins), utheta(nbins,ndim))

		allocate(prestheta_local(nbins,ndim), visctheta_local(nbins,ndim))
		allocate(prestheta(nbins,ndim), visctheta(nbins,ndim))

		theta_bin = zero
		ptheta_local = zero
		utheta_local = zero
		prestheta_local = zero
		visctheta_local = zero
		theta_count_local = 0

		ptheta = zero
		utheta = zero
		prestheta = zero
		visctheta = zero
		theta_count = 0

		call compute_omega

		do itheta=1, nbins
			theta_bin(itheta) = (itheta-0.5) * dtheta
		enddo


		do m=1, nbody
			if (myid_particles(m)) then

				da=4.*pi*(radbdy(m)*dx)**2./real(nbnd,prcn)
     
				iphs = 1
				nbnd = phase_array(iphs)%nbnd
				nullify(bndarray)
				bndarray => phase_array(iphs)%bndpts


				do l=1, nbnd
					if (.not.part_array(m)%if_drag(L)) goto 100
					rad = zero
					do n=1, ndim
						xl(n)=xc(m,n)+ bndarray(n,l)*radbdy(m)

						ul(n)=zero
						ppll(n)=zero
				 	enddo

					rad  = DSQRT(dot_product(bndarray(1:ndim,l),bndarray(1:ndim,l)))
					rad2 = DSQRT(dot_product(bndarray(2:ndim,l),bndarray(2:ndim,l)))-small_number

					do n=1, ndim
						vcellb(n) = floor(xl(n))
						pcellb(n) = floor(xl(n))
					enddo

					i_have_the_point = point_in_this_domain(vcellb, index_out, xl, position_out)

					if (i_have_the_point) then
						vcellb(:) = index_out(:)
						xl(:) = position_out(:)

						pcellb(:) = vcellb(:)

						pl=zero
						call interpolate_pdata(pcellb, xl, ppll, pl, l)

						pl = pl + (xl(1) - dble(global_n(1))/2)/global_n(1) * doml(1) * mpg(1)


						call interpolate_udata(vcellb, xl, ib, ie, jb, je, kb, ke, ul, nll, onll, dfll, 0, m, l, onew) 

						vort(:) = zero 
						do k = 1, onew 
							do j = 1, onew
								do i = 1, onew
									ii = ib+i-1
									jj = jb+j-1
									kk = kb+k-1
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

						xtemp = -bndarray(1,l)
						theta = acos(abs(xtemp)/rad)
						if (xtemp<zero) theta = pi - theta

						!rad2 = rad*sin(theta)
						if (rad2>small_number) then
							ytemp = bndarray(2,l)
							ztemp = bndarray(3,l)

							phi = acos(abs(ytemp)/rad2)
							if (ytemp>zero .and. ztemp>zero) then
							elseif (ytemp<zero .and. ztemp>zero) then
								phi = pi-phi
							elseif (ytemp<zero .and. ztemp<zero) then
								phi = pi+phi
							elseif (ytemp>zero .and. ztemp<zero) then
								phi = 2*pi-phi
							endif
						else
							phi = 0
						endif

						itheta = floor(theta / dtheta) + 1
						if (itheta > nbins) itheta = nbins

						theta_count_local(itheta) = theta_count_local(itheta) + 1
						ptheta_local(itheta) = ptheta_local(itheta) + pl

						call vec_cart_to_sphr(ul,theta,phi)
						utheta_local(itheta, :) = utheta_local(itheta, :) + ul(:)


						!---------------------------------------------------------------
						!     calculate the viscous and pressure components separately
						!---------------------------------------------------------------
						do n=1, ndim
							prestheta_local(itheta,n) = prestheta_local(itheta,n) - pl * cd(n,l) * da

							visctheta_local(itheta,n) = visctheta_local(itheta,n) + df(l,n) * vis * da
						enddo
					endif
100				continue
				enddo
			endif
		enddo

		GLOBAL_DOUBLE_SUM(ptheta_local, ptheta, nbins, comm_cart_2d)
		GLOBAL_DOUBLE_SUM(utheta_local, utheta, nbins*ndim, comm_cart_2d)

		GLOBAL_DOUBLE_SUM(prestheta_local, prestheta, nbins*ndim, comm_cart_2d)
		GLOBAL_DOUBLE_SUM(visctheta_local, visctheta, nbins*ndim, comm_cart_2d)

		GLOBAL_LONGINT_SUM(theta_count_local, theta_count, nbins, comm_cart_2d)

		do itheta=1, nbins
			if (theta_count(itheta)>0) then
				ptheta(itheta) = ptheta(itheta) / theta_count(itheta)
				utheta(itheta, :) = utheta(itheta, :) / theta_count(itheta)
			endif
		enddo

		ptheta(:) = ptheta(:) / (vis*umeanslip/dia_phys)
		utheta(:, :) = utheta(:, :) / umeanslip

		prestheta(:,:) = prestheta(:,:) / (.5 * rhof * umeanslip**2 * pi * dia_phys**2/4)  ! (3*pi*dia_phys*vis*(1-maxvolfrac)*umeanslip)
		visctheta(:,:) = visctheta(:,:) / (.5 * rhof * umeanslip**2 * pi * dia_phys**2/4) !(3*pi*dia_phys*vis*(1-maxvolfrac)*umeanslip)

		if (I_AM_NODE_ZERO) then
			filename = trim(run_name)//"_p_u_theta"
			unitnum = 1
			call instant_file_opener(filename,unitnum,.true.)

			do itheta=1, nbins			
				if (theta_count(itheta)>0) write (unitnum, "(1d15.7,1i,10d15.7)") theta_bin(itheta), int(theta_count(itheta)), ptheta(itheta), utheta(itheta,:), prestheta(itheta,:), visctheta(itheta,:)
			enddo
			close(unitnum)
			write (*,*) "PU_vs_THERA IS DONE..."
			call screen_separator(80,'t')
		endif
	end subroutine pressure_velocity_vs_theta


	subroutine vec_cart_to_sphr(vec,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: vec(ndim)

		real(8) :: tmp(3)
		real(8) :: rotate(3,3)

		integer :: dim1,dim2,dim3,dim4

		!ROTATION OF THE STRESS TENSOR ALIGNED WITH THE SPHERICAL COORDINATE
		rotate(1,1) =-cos(th)
		rotate(2,1) = sin(th) * cos(phi)
		rotate(3,1) = sin(th) * sin(phi)

		rotate(1,2) = sin(th)
		rotate(2,2) = cos(th) * cos(phi)
		rotate(3,2) = cos(th) * sin(phi)

		rotate(1,3) = 0d0
		rotate(2,3) =-sin(th) * sin(phi)
		rotate(3,3) = sin(th) * cos(phi)

		tmp = 0d0
		do dim1=1, ndim
			do dim2=1, ndim
				tmp(dim1) = tmp(dim1) + rotate(dim2,dim1) * vec(dim2)
			enddo
		enddo
		vec(:) = tmp(:)
	end subroutine vec_cart_to_sphr

	subroutine vec_sphr_to_cat(vec,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: vec(ndim)

		real(8) :: tmp(3)
		real(8) :: rotate(3,3)

		integer :: dim1,dim2,dim3,dim4

		!ROTATION OF THE STRESS TENSOR ALIGNED WITH THE SPHERICAL COORDINATE
		rotate(1,1) =-cos(th)
		rotate(2,1) = sin(th) * cos(phi)
		rotate(3,1) = sin(th) * sin(phi)

		rotate(1,2) = sin(th)
		rotate(2,2) = cos(th) * cos(phi)
		rotate(3,2) = cos(th) * sin(phi)

		rotate(1,3) = 0d0
		rotate(2,3) =-sin(th) * sin(phi)
		rotate(3,3) = sin(th) * cos(phi)

		tmp = 0d0
		do dim1=1, ndim
			do dim2=1, ndim
				tmp(dim1) = tmp(dim1) + rotate(dim1,dim2) * vec(dim2)
			enddo
		enddo
		vec(:) = tmp(:)
	end subroutine vec_sphr_to_cat


end module mypost_process

