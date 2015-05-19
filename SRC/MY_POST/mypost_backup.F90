MOdule mypost_process
#include "../FLO/ibm.h"
	use precision
	use global_data
	use nlmainarrays , only : ubcp, pbcp
	use fftw_interface
	use constants
	use nlarrays
	use randomno
	use general_funcs
	use string_funcs
	use boundary_condition
	use dependent_functions
	use bcsetarrays
	use epsarray
	use postproc_funcs
	use clustermod
!	use string_funcs
!	use init_turb
	implicit none
!	private
!	public :: dissip_continuous, compute_sijsij

	logical, allocatable, dimension(:,:,:) :: far_field, fluid_atijkp
	logical :: near_particle_checked=.false.
	real(prcn) :: near_particle_rad = 1.05

	real(prcn) :: tke, ufmean_mag, usmean1_mag, usmean2_mag, mix_mean_vel_mag, mix_mean_slip_mag, sijsij, sijsij_pure, tke_pure, source, dissipation, int_tke_pure

	real(prcn), allocatable :: mix_mean_slip(:), mean_spec_vel_mag(:), &
		& mix_mean_vel(:), grant(:), tke_s_pure(:), acc_mean(:), vel_mean(:)

!	integer :: nvar, nvar1, nvar2, 
	integer :: line=1

!	type :: slice_type
!		integer :: num, dy
!		integer, allocatable :: ind(:), fluidnum(:), solidnum(:)
!		real(prcn), allocatable :: ufmean(:,:), mpg(:,:)
!		real(prcn), allocatable :: areaf(:)
!	end type slice_type
!	type(slice_type),allocatable:: slice(:)
contains
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

			mix_mean_vel(:) =  mix_mean_vel(:) + phase_array(iphs)%mean_spec_vel(:) * phase_array(iphs)%volfracg / maxvolfrac


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

	subroutine rotate_vec_forward(direct_cos, vec_in, vec_out, n)
		implicit none
		integer, intent(in) :: n
		real(prcn), intent(in) :: direct_cos(n,n), vec_in(n)
		real(prcn), intent(out) :: vec_out(n)
		integer :: dim1, dim2

		vec_out = zero

		do dim1=1, n
			do dim2=1, n
				vec_out(dim1) = vec_out(dim1) + direct_cos(dim2,dim1) * vec_in(dim2)
			enddo
		enddo
	end subroutine rotate_vec_forward

	subroutine rotate_vec_backward(direct_cos, vec_in, vec_out, n)
		implicit none
		integer, intent(in) :: n
		real(prcn), intent(in) :: direct_cos(n,n), vec_in(n)
		real(prcn), intent(out) :: vec_out(n)
		integer :: dim1, dim2

		vec_out = zero

		do dim1=1, n
			do dim2=1, n
				vec_out(dim1) = vec_out(dim1) + direct_cos(dim1,dim2) * vec_in(dim2)
			enddo
		enddo
	end subroutine rotate_vec_backward

	subroutine post_part_stat
		implicit none

		type :: slice_type
			integer :: num
			real(prcn) :: start, end
			real(prcn) :: usmean(ndim), ufmean(ndim), rij(ndim, ndim), bij(ndim, ndim), grant, dens
			integer, allocatable :: list(:)
			real(prcn), allocatable :: gofr(:)
		end type slice_type
		type(slice_type),allocatable:: slice_y(:), slice_y_tmp(:)
		integer :: slice_n

		real(prcn) :: max_overlap, grant
		integer :: ibody, ibin, count, io, islice
		real(prcn), allocatable :: xc_2d_tmp(:,:)
		logical, allocatable :: contact(:,:)
		integer, allocatable :: plane_color(:)
		character*50 filename1, filename2, filename3, filename4
		logical :: filexist

		first_pass = .true.
		if (I_AM_NODE_ZERO) then
			write (*,*) "IN POST_PART_STAT"

			if (.not.allocated(plane_color)) allocate(plane_color(nbody))
			plane_color = 0

			if (allocated(rad_bin)) deallocate(rad_bin)
			allocate(rad_bin(nbins))
			if (allocated(gofr_avg)) deallocate(gofr_avg)
			allocate(gofr_avg(nbins))

			allocate(contact(nbody,nbody))

!			filename2 = TRIM(RUN_NAME)//"_part_position.dat"
!			open (unit=10002, file=trim(filename2), status="replace", action="write")

!			filename3 = TRIM(RUN_NAME)//'_gofr.dat'
!			open (unit=10003, file=trim(filename3), status="replace", action="write")

			filename4 = TRIM(RUN_NAME)//"_cluster"

			filename1 = TRIM(RUN_NAME)//"_part_info.rst"
			inquire (file=trim(filename1), exist=filexist)
			if (filexist) then
				write (*,*) "FILE "//filename1//" EXISTS"
				open (unit=10001, file=trim(filename1), status="old", action="read", form="unformatted")
			else
				write (*,*) "FILE  "//filename1//" DOES NOT EXIST"
				goto 100
!				stop
			endif

			count=0
			do
				count = count+1

				read(10001,iostat=io) t
				read(10001,iostat=io) xc(1:nbody, 1:ndim)
				read(10001,iostat=io) velbdy(1:nbody, 1:ndim)
				read(10001,iostat=io) force(1:nbody, 1:ndim)
				read(10001,iostat=io) pres(1:nbody, 1:ndim)
				read(10001,iostat=io) visc(1:nbody, 1:ndim)
				read(10001,iostat=io) contact_force(1:nbody, 1:ndim)
				read(10001,iostat=io) frame_vel(1:ndim)
				read(10001,iostat=io) frame_accln(1:ndim)
				read(10001,iostat=io) ufmean(1:ndim)

				if (io<0) then
					write (*,*) "END OF FILE ", trim(filename1)

#if 0
					open (unit=1, file="part_vel.dat", status="replace", action="write")
					write (1,"(5d15.7)") (xc(ibody,:), radbdy(ibody), sqrt(dot_product(velbdy(ibody,:),velbdy(ibody,:))), ibody=1,nbody)
					close(1)
#endif
					exit

				elseif(io==0) then
					if (mod(count,10)==1) then
						write (*,"(1a,1i,1a)") "IN THE ", count, "th ZONE"

						xc(:,:) = xc(:,:)/dbydx
						if (first_pass) radbdy(:) = radbdy(:)/dbydx

						mx1 = lybyd !mx1/dbydx
						nx = lybyd !nx/dbydx
						my = lybyd !my/dbydx
						mz = lybyd !mz/dbydx

!						if (first_pass) call calc_plane_color(nbody, xc, plane_color)
						call total_grant
						call form_slice
						call compute_2d_values
!						call make_output


!						call calculate_gofr_homog(nbody, xc, contact, my, mx+1, nbins, .true., gofr_avg, rad_bin, max_overlap)
!						call find_clusters(xc, contact, nbody, lybyd, filename4, .true.)

!						write (10002,*) "zone"
!						write (10002,"(4d15.7,2i)") ((xc(ibody,:), radbdy(ibody), aggregate_num, plane_color(ibody)), ibody=1, nbody)
!						write (10002,"(4d15.7,1i)") ((xc(ibody,:), radbdy(ibody), aggregate_num(ibody)), ibody=1, nbody)

!						write (10003,*) "zone"
!						write (10003,"(2d15.7)") ((rad_bin(ibin), gofr_avg(ibin)), ibin=1,nbins)

#if 0
						do islice=1, slice_n
							if (slice_y(islice)%num>0) then
								deallocate(slice_y(islice)%list)
								deallocate(slice_y(islice)%gofr)
								slice_y(islice)%num = 0
								slice_y(islice)%usmean = zero
								slice_y(islice)%rij = zero
								slice_y(islice)%bij = zero
								slice_y(islice)%grant = zero
							endif
						enddo
#endif
						first_pass = .false.
					endif
				endif

!if (count==3000) exit
			enddo			



100		continue
#if 0
			if (.not.filexist) then			
				call calc_plane_color(nbody, xc, plane_color)
				call calculate_gofr_homog(nbody, xc, contact, my, mx+1, nbins, .true., gofr_avg, rad_bin, max_overlap)

				write (10002,*) "zone"
				write (10002,"(4d15.7,1i)") ((xc(ibody,:), radbdy(ibody), plane_color(ibody)), ibody=1, nbody)

				write (10003,*) "zone"
				write (10003,"(2d15.7)") ((rad_bin(ibin), gofr_avg(ibin)), ibin=1,nbins)
			endif
#endif
			write (*,*) "OUT OF POST_PART_STAT"
!			deallocate(plane_color, contact)

			close(1001)
!			close(1002)
!			close(1003)
		endif

	contains

		subroutine form_slice
			implicit none
			integer :: ibody, int_tmp, islice

			if (.not.allocated(slice_y)) then
				if ( mod(my,slice_dx*2)/=0) then
					write (*,*) "CHOOSE ANOTHER VALUE FOR SLICE_DX. THERE IS NO STRAIGHTFORWARD LAYER FORMATION"
					stop
				endif

				slice_n = my / slice_dx/2

				allocate(slice_y(slice_n), slice_y_tmp(slice_n))
				do islice=1, slice_n
					slice_y(islice)%start = (islice-1) * slice_n
					slice_y(islice)%end   =  islice    * slice_n

					slice_y_tmp(islice)%start = (islice-1) * slice_n
					slice_y_tmp(islice)%end   =  islice    * slice_n
				enddo

				do islice=1, slice_n
					allocate(slice_y_tmp(islice)%list(nbody))
				enddo
			endif

			do islice=1, slice_n
				slice_y_tmp(islice)%num = 0
				slice_y_tmp(islice)%list(:) = 0
			enddo

			do ibody=1, nbody
				islice = xc(ibody,2) / slice_dx + 1
				if (islice>2*slice_n) islice=2*slice_n

				if (islice>slice_n) islice = 2*slice_n - islice + 1


!write (*,*) ibody, islice
!if (ibody==25) write (*,"(3d15.7)") xc(ibody,:)

				if (islice<1.or.islice>slice_n) then
					write (*,*) "ISLICE IS OUT OF RANGE, ", islice, ibody

					if (islice<1) isclice = 1
					if (islice>slice_n) isclice = slice_n
!write (*,"(3d15.7)") xc(ibody,:)
				else
					slice_y_tmp(islice)%num = slice_y_tmp(islice)%num+1
					slice_y_tmp(islice)%list(slice_y_tmp(islice)%num) = ibody
				endif
			enddo

			do islice=1, slice_n
!				if (allocated(slice_y(islice)%list)) deallocate(slice_y(islice)%list)
				int_tmp = slice_y_tmp(islice)%num
				slice_y(islice)%num = int_tmp
				if (int_tmp>0) then
					allocate(slice_y(islice)%list(int_tmp))
					slice_y(islice)%list(1:int_tmp) = slice_y_tmp(islice)%list(1:int_tmp)

					allocate(slice_y(islice)%gofr(nbins))
					slice_y(islice)%gofr(:) = zero

					slice_y(islice)%usmean = zero
					slice_y(islice)%ufmean = zero
					slice_y(islice)%rij = zero
					slice_y(islice)%bij = zero
					slice_y(islice)%grant = zero
					slice_y(islice)%dens = zero
				endif
			enddo

			do islice=1, slice_n
				deallocate(slice_y_tmp(islice)%list)
			enddo
			deallocate(slice_y_tmp)

		end subroutine form_slice

		subroutine compute_2d_values
			implicit none

			integer :: ibody, islice, i, j, k, fcount
			real(prcn), allocatable :: xc_2d_tmp(:,:), rad_tmp(:)
			real(prcn) :: norm_dens, vec1(ndim)

			if (allocated(rad_bin)) deallocate(rad_bin)
			allocate(rad_bin(nbins))

			norm_dens = nbody / (lybyd*lybyd*lybyd/2)

			do islice=1, slice_n
				if (slice_y(islice)%num>0) then
					allocate(xc_2d_tmp(slice_y(islice)%num,2))
					allocate(rad_tmp(slice_y(islice)%num))
					do ibody=1, slice_y(islice)%num
						xc_2d_tmp(ibody,1) = xc(slice_y(islice)%list(ibody),1)
						xc_2d_tmp(ibody,2) = xc(slice_y(islice)%list(ibody),3)
						rad_tmp(ibody) = radbdy(slice_y(islice)%list(ibody))
					enddo

					slice_y(islice)%dens = slice_y(islice)%num / (lybyd*lybyd*lybyd/slice_n) / norm_dens


!					call calc_gofr_2d(slice_y(islice)%num, xc_2d_tmp, int(lybyd), int(lybyd), nbins, rad_tmp, rad_bin, slice_y(islice)%gofr)
					call garn_tmp(slice_y(islice)%num, slice_y(islice)%list, slice_y(islice)%usmean, &
					&             slice_y(islice)%rij, slice_y(islice)%bij, slice_y(islice)%grant)

!					vec1(:) = zero
!					fcount = 0
!					do i=slice_y(islice)%start, slice_y(islice)%end-1
!						do j=1, my
!							do k=1, mz
!								if (fluid_atijk(i,j,k)) then
!									vec1(:) = vec1(:)+ubcp(:)
!									fcount = fcount+1
!								endif
!							enddo
!						enddo
!					enddo
!					slice_y(islice)%ufmean(:) = vec1(:) / fcount

					deallocate(xc_2d_tmp, rad_tmp)
				endif
			enddo
		end subroutine compute_2d_values

		subroutine garn_tmp(n, list, means, rij, bij, grant)
			implicit none
			integer, intent(in) :: n
			integer, intent(in) :: list(n)
			real(prcn), intent(out) :: means(ndim), rij(ndim,ndim), bij(ndim,ndim), grant
			integer :: ibody, dim1, dim2
			real(prcn) :: xi, eta

			means = zero
			grant = zero

			do ibody=1, n
				means(:) = means(:) + velbdy(list(ibody),:)
			enddo
			means(:) = means(:)/n

			do ibody=1, n
				do dim1=1, ndim
					do dim2=1, ndim
						rij(dim1,dim2) = rij(dim1,dim2) + (velbdy(list(ibody),dim1)-means(dim1)) * (velbdy(list(ibody),dim2)-means(dim2))
					enddo
				enddo
			enddo
			rij = rij/n
			grant = (rij(1,1) + rij(2,2) + rij(3,3)) / 3

			call calc_anisotropy(rij, bij, xi, eta)
		end subroutine garn_tmp

		subroutine total_grant
			implicit none

			integer :: ibody, idim

			usmean(:) = zero
			grant = zero
			do ibody = 1, nbody
				usmean(:) = usmean(:) + velbdy(ibody,:)
			end do
			usmean(:) = usmean(:)/nbody

			do ibody=1, nbody
				do idim = 1, ndim
					grant = grant + (velbdy(ibody,idim)-usmean(idim))**2.d0
				end do
			end do
			grant = grant/(three*nbody)
		end subroutine total_grant


		subroutine make_output
			implicit none
			integer :: islice, ibin, i, j
			character*50 filename1, filename2, filename3, filename4, filename5
			real(prcn) :: tmp1, tmp2, tmp3, deltan, deltau, eta, ybar, nc, ns, fdens, uc, us, fu, xs, xe, ys, ye, f, vec1(ndim)
			integer :: iln, ihn, ilu, ihu

#if 0
			filename1=trim(run_name)//"_gofr_2d.dat"
			if (first_pass) then
				open (unit=2, file=trim(filename1), status="replace", action="write")
			else
				open (unit=2, file=trim(filename1), status="old", action="write", position="append")
			endif
!			write (2,"(1a,1i,1a,1i,1a)") "zone i= ", nbins, " j= ", slice_n, " F=point"

			do islice=1, slice_n
				write (2,*) 'ZONE T = "', t/t_conv/(1-maxvolfrac), '",'
				if (slice_y(islice)%num>0) then
					do ibin=1, nbins
						write (2,"(2i,3d15.7)") ibin, islice, rad_bin(ibin), (islice-0.5) * my/2/slice_n, slice_y(islice)%gofr(ibin)
					enddo
				else
					do ibin=1, nbins
						write (2,"(2i,3d15.7)") ibin, islice, rad_bin(ibin), (islice-0.5) * my/2/slice_n, zero
					enddo
				endif
			enddo
			close(2)
#endif

			filename2=trim(run_name)//"_slice_values.dat"
			if (first_pass) then
				open (unit=3, file=trim(filename2), status="replace", action="write")
			else
				open (unit=3, file=trim(filename2), status="old", action="write", position="append")
			endif
			write (3,*) 'ZONE T = "', t/t_conv/(1-maxvolfrac), '",'

			do islice=1, slice_n
!				if (slice_y(islice)%num>0) then
					tmp1 = sqrt(dot_product(slice_y(islice)%usmean(:),slice_y(islice)%usmean(:)))
					tmp2 = ucharmod/(1-maxvolfrac)   !sqrt(dot_product(ufmean(:),ufmean(:)))
!					vec1(:) = slice_y(islice)%ufmean(:)-slice_y(islice)%usmean(:)
!					tmp3 = sqrt(dot_product(vec1,vec1))
!					slice_y(islice)%re = slice_y(islice)%num / ((slice_y(islice)%end-slice_y(islice)%start)*my*mz*dx**3) * pi*dia_phys**3/6 * tmp3*dia_phys/vis 

!					write (3,"(39d15.7)") (islice-0.5) * my/2/slice_n, slice_y(islice)%dens, slice_y(islice)%usmean(:)/tmp2, slice_y(islice)%usmean(:)/tmp3, tmp1/tmp2, tmp1/tmp3, tmp2/tmp3, slice_y(islice)%rij(:,:)/grant, ((slice_y(islice)%rij(i,j)/slice_y(islice)%grant, j=1,ndim), i=1,ndim), slice_y(islice)%grant/grant, slice_y(islice)%bij(:,:)

					write (3,"(34d15.7)") (islice-0.5) * my/2/slice_n, slice_y(islice)%dens, slice_y(islice)%usmean(:)/tmp2, tmp1/tmp2, slice_y(islice)%rij(:,:)/grant, ((slice_y(islice)%rij(i,j)/slice_y(islice)%grant, j=1,ndim), i=1,ndim), slice_y(islice)%grant/grant, slice_y(islice)%bij(:,:)

!				endif
			enddo
			close(3)


			filename3=trim(run_name)//"_slice_values_self_similar_n.dat"
			filename4=trim(run_name)//"_slice_values_self_similar_u.dat"
			filename5=trim(run_name)//"_slice_values_delta.dat"
			if (first_pass) then
				open (unit=3, file=trim(filename3), status="replace", action="write")
				open (unit=4, file=trim(filename4), status="replace", action="write")
				open (unit=5, file=trim(filename5), status="replace", action="write")
			else
				open (unit=3, file=trim(filename3), status="old", action="write", position="append")
				open (unit=4, file=trim(filename4), status="old", action="write", position="append")
				open (unit=5, file=trim(filename5), status="old", action="write", position="append")
			endif

			iln = 0
			ihn = 0
			ilu = 0
			ihu = 0

			do islice=1, slice_n-1
				tmp1 = abs(slice_y(islice)%dens-minval( slice_y(:)%dens )) / &
					&	abs(maxval( slice_y(:)%dens )-minval( slice_y(:)%dens ))

				tmp2 = abs(slice_y(islice+1)%dens-minval( slice_y(:)%dens )) / &
					&	abs(maxval( slice_y(:)%dens )-minval( slice_y(:)%dens ))

				if (tmp1>=0.9 .and. 0.9>=tmp2) iln = islice
				if (tmp1>=0.1 .and. 0.1>=tmp2) ihn = islice

				tmp1 = abs(slice_y(islice)%usmean(1)-minval( slice_y(:)%usmean(1) )) / &
					&	abs(maxval( slice_y(:)%usmean(1) )-minval( slice_y(:)%usmean(1) ))

				tmp2 = abs(slice_y(islice+1)%usmean(1)-minval( slice_y(:)%usmean(1) )) / &
					&	abs(maxval( slice_y(:)%usmean(1) )-minval( slice_y(:)%usmean(1) ))


				if (tmp1<=0.1 .and. 0.1<=tmp2) ilu = islice
				if (tmp1<=0.9 .and. 0.9<=tmp2) ihu = islice
			enddo

!			if (count>=100) then
				if (iln>0.and.ihn>0.and.iln/=ihn) then
					write (3,*) 'ZONE T = "', t/t_conv/(1-maxvolfrac), '",'

					ys = 0.9 * abs(maxval( slice_y(:)%dens )-minval( slice_y(:)%dens )) + minval( slice_y(:)%dens )
					ye = 0.1 * abs(maxval( slice_y(:)%dens )-minval( slice_y(:)%dens )) + minval( slice_y(:)%dens )

					xs = iln-0.5 + (ys-slice_y(iln)%dens) * (iln+1-iln) / (slice_y(iln+1)%dens-slice_y(iln)%dens)
					xe = ihn-0.5 + (ye-slice_y(ihn)%dens) * (ihn+1-ihn) / (slice_y(ihn+1)%dens-slice_y(ihn)%dens)

					nc = 0.5*(ye+ys)
					ns = (ys-ye)

					do islice=iln, ihn+1
						deltan = (xe-xs)
						ybar = 0.5*(xs+xe)

						if (islice==iln) then
							eta = (xs - ybar) / deltan
							f = ys
						elseif (islice==ihn+1) then
							eta = (xe - ybar) / deltan
							f = ye
						else
							eta = ((islice-0.5) - ybar) / deltan
							f = slice_y(islice)%dens
						endif

						fdens = (f - nc)/ns
				
						write (3,"(2d15.7)") eta, fdens
					enddo
				endif

				if (ilu>0.and.ihu>0.and.ilu/=ihu) then
					write (4,*) 'ZONE T = "', t/t_conv/(1-maxvolfrac), '",'

					ys = 0.1 * abs(maxval( slice_y(:)%usmean(1) )-minval( slice_y(:)%usmean(1) )) + minval( slice_y(:)%usmean(1) )
					ye = 0.9 * abs(maxval( slice_y(:)%usmean(1) )-minval( slice_y(:)%usmean(1) )) + minval( slice_y(:)%usmean(1) )

					xs = ilu-0.5 + (ys-slice_y(ilu)%usmean(1)) * (ilu+1-ilu) / (slice_y(ilu+1)%usmean(1)-slice_y(ilu)%usmean(1))
					xe = ihu-0.5 + (ye-slice_y(ihu)%usmean(1)) * (ihu+1-ihu) / (slice_y(ihu+1)%usmean(1)-slice_y(ihu)%usmean(1))

					nc = 0.5*(ye+ys)
					ns = (ye-ys)

					do islice=ilu, ihu+1
						deltau = (xe-xs)
						ybar = 0.5*(xs+xe)

						if (islice==ilu) then
							eta = (xs - ybar) / deltau
							f = ys
						elseif (islice==ihu+1) then
							eta = (xe - ybar) / deltau
							f = ye
						else
							eta = ((islice-0.5) - ybar) / deltau
							f = slice_y(islice)%usmean(1)
						endif

						fu = (f - nc)/ns
				
						write (4,"(2d15.7)") eta, fu
					enddo
				endif

				write (5,"(4d15.7)") t/t_conv/(1-maxvolfrac), deltan, sqrt(grant), deltan/sqrt(grant)
!			endif
			close(3)
			close(4)
			close(5)
		end subroutine make_output


	end subroutine post_part_stat

	subroutine calc_plane_color(nbody, xc, color)
		implicit none
		integer, intent(in) :: nbody
		real(prcn), intent(in) :: xc(nbody,ndim)
		integer, intent(out) :: color(:)
		integer :: ibody, length

		length = size(color,1)
		if (length/=nbody) then
			write (*,*) "THE LENGTH DOES NOT MATCH NBODY IN CALC_PLANE_COLOR, CHECK IT..."
			stop
		endif

		do ibody=1, nbody
			if (xc(ibody,2)/dbydx<=dia_phys) then
				color(ibody) = 2
			elseif (lybyd/4-dia_phys<=xc(ibody,2)/dbydx .and. xc(ibody,2)/dbydx<=lybyd/4) then
				color(ibody) = 3
			elseif (lybyd*3/4<=xc(ibody,2)/dbydx .and. xc(ibody,2)/dbydx<=lybyd*3/4+dia_phys) then
				color(ibody) = 4
			elseif (lybyd-dia_phys<=xc(ibody,2)/dbydx .and. xc(ibody,2)/dbydx<=lybyd) then
				color(ibody) = 5
			else
				color(ibody) = 1
			endif
		enddo
	end subroutine calc_plane_color


	subroutine relative_acceleration
		implicit none
		integer :: i, j, ibin, idim, nvar, nvar1, nvar2
		real(prcn) :: rmax, dr, r, rij(ndim), rij_unit(ndim), acc_dif(ndim), acc_dif_mag, mass, vec1(ndim), sep, norm_factor, tmp
		real(prcn), allocatable :: acc_hydro_rad(:), acc_hydro_vec(:,:) !, acc_part(:,:), acc_part_rad(:), acc_part_vec(:,:)
		integer(8), allocatable :: num_bin(:)
		real(prcn), allocatable :: acc_fluc(:,:), vel_fluc(:,:), acc_var(:,:), vel_var(:,:), acc_fluc_meanf(:,:), acc_var_meanf(:,:)
		character*80 filename1
		logical :: filexist

		if (I_AM_NODE_ZERO) then
			write (*,*) "IN RELATIVE_ACCELERATION"
			if (.not. post_no_flow_mem_alloc) then
				call compute_mix_mean_vel

				allocate(acc_fluc(nbody,ndim), vel_fluc(nbody,ndim), acc_var(nphases,ndim), vel_var(nphases,ndim))
				allocate(acc_fluc_meanf(nbody,ndim), acc_var_meanf(nphases,ndim))
!				call AiVi(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)

				if (allocated(acc_hydro_rad)) deallocate(acc_hydro_rad)
				allocate(acc_hydro_rad(nbins), acc_hydro_vec(nbins,ndim))
!				allocate(acc_part(nbins,ndim), acc_part_vec(nbins,ndim), acc_part_rad(nbins))

				if (allocated(rad_bin)) deallocate(rad_bin)
				allocate(rad_bin(nbins), num_bin(nbins))

				acc_hydro_rad = zero
				acc_hydro_vec = zero

!				acc_part = zero
!				acc_part_rad = zero
!				acc_part_vec = zero

				rad_bin = zero
				num_bin = 0

				rmax  = lybyd*dbydx / 2
				dr = rmax / nbins
				do ibin=1, nbins
					rad_bin(ibin) = (ibin-0.5) * dr
				enddo

#if 0
				do i=1, nbody-1
					do j=i+1, nbody
						rij(:) = xc(i,:) - xc(j,:)
						do idim=1, ndim
							if (rij(idim) > rmax) then
								rij(idim) = rij(idim) - 2*rmax
							elseif (rij(idim) < -rmax) then
								rij(idim) = rij(idim) + 2*rmax
							endif
						enddo
						r = sqrt(dot_product(rij,rij))
						rij_unit(:) = rij(:) / r

						sep = r/dbydx - dia_phys
						if (sep<d0_rat*dia_phys) sep = d0_rat*dia_phys
						acc_part(i,:) = acc_part(i,:) + hamaker*dia_phys / (12*sep**2) * rij_unit(:) / (nbody-1)
						acc_part(j,:) = acc_part(j,:) - hamaker*dia_phys / (12*sep**2) * rij_unit(:) / (nbody-1)

if (i==157.or.j==157) then
!write (*,"(2d15.7)") r, sep
!write (*,"(3d15.7)") rij_unitl
!write (*,"(4d15.7)") hamaker, dia_phys, (12*sep**2), nbody-1
!write (*,"(3d15.7)") dia_phys / (12*sep**2) / nbody-1
!read (*,*)
write (*,"(1i, 3d15.7)") i, acc_part(i,:)
endif

					enddo
				enddo
#endif
				mass = rhos*pi*dia_phys**3.d0/6.d0
				do i=1, nbody-1
					do j=i+1, nbody
						rij(:) = xc(i,:) - xc(j,:)
						do idim=1, ndim
							if (rij(idim) > rmax) then
								rij(idim) = rij(idim) - 2*rmax
							elseif (rij(idim) < -rmax) then
								rij(idim) = rij(idim) + 2*rmax
							endif
						enddo
						r = sqrt(dot_product(rij,rij))

						if (r<=rmax) then
							ibin = int(r/rmax * nbins) +1
							if (ibin>nbins) ibin = nbins
							if (rad_bin(ibin)/dbydx < dia_phys) ibin = ibin+1

							rij_unit(:) = rij(:) / r
							acc_dif(:) = (force(i,:) - force(j,:))
							acc_dif_mag = dot_product(acc_dif(:), rij_unit(:))

!if (i==167.and.j==168) then
!write (*,*) ibin, acc_dif_mag
!write (*,"(9d15.7)") acc_dif(:), force(i,:), force(j,:)
!write (*,"(4d15.7)") rij(:), r
!endif

							acc_hydro_vec(ibin,:) = acc_hydro_vec(ibin,:) + acc_dif(:)
							acc_hydro_rad(ibin)   = acc_hydro_rad(ibin)   + acc_dif_mag

!							acc_dif(:) = (acc_part(i,:) - acc_part(j,:))
!							acc_dif_mag = dot_product(acc_dif(:), rij_unit(:))

!							acc_part_vec(ibin,:) = acc_part_vec(ibin,:) + acc_part(i,:)-acc_part(j,:)
!							acc_part_rad(ibin)   = acc_part_rad(ibin)   + acc_dif_mag

							num_bin(ibin) = num_bin(ibin) + 1
						endif
					enddo
				enddo

				do ibin=1, nbins
					if (num_bin(ibin)>0) then
						acc_hydro_vec(ibin,:) = acc_hydro_vec(ibin,:) / num_bin(ibin)
						acc_hydro_rad(ibin)   = acc_hydro_rad(ibin)   / num_bin(ibin)

!						acc_part_vec(ibin,:) = acc_part_vec(ibin,:) / num_bin(ibin)
!						acc_part_rad(ibin)   = acc_part_rad(ibin)   / num_bin(ibin)
					endif
				enddo

				norm_factor = 3.d0*pi*vis*mix_mean_slip_mag*(1-maxvolfrac)*dia_phys


				filename1=trim(run_name)//"_relative_acceleration"
				if (from_post) then
					filename1=trim(filename1)//"_post.dat"
					open (unit=1, file=trim(filename1), status="replace", action="write")
				else
					filename1=trim(filename1)//".dat"
					if (irestart==0.and.first_pass) then
						open (unit=1, file=trim(filename1), status="replace", action="write")
					else
						inquire (file=trim(filename1), exist=filexist)
						if (.not.filexist) then
							open (unit=1, file=trim(filename1), status="replace", action="write")
						else
							open (unit=1, file=trim(filename1), status="old", action="write", position="append")
						endif
					endif
				endif

				write (1,*) "zone"
				do ibin=1, nbins
					tmp = sqrt(acc_hydro_vec(ibin,1)**2+acc_hydro_vec(ibin,2)**2+acc_hydro_vec(ibin,3)**2)

					write (1,"(1f10.4,8d15.7)") rad_bin(ibin)/dbydx, acc_hydro_rad(ibin)/norm_factor, acc_hydro_vec(ibin,:)/norm_factor, tmp/norm_factor
!, acc_part_rad(ibin)/norm_factor, acc_part_vec(ibin,:)/norm_factor
				enddo
				close(1)

!				deallocate(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)
!				deallocate(acc_hydro_rad, acc_hydro_vec, acc_part_rad, acc_part_vec, rad_bin)
			else
				nvar  = 5
				nvar1 = 1
				nvar2 = 4
				line = nbins

				filename1 = "_relative_acceleration.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line)
			endif
		endif
	end subroutine relative_acceleration


	SUBROUTINE flow_snapshot2
		Use nlarrays , only : ur1, uf1
		use nlmainarrays, only : ubcp, pbcp
		USE dem_mod, only : is_mobile, des_pos_new, des_radius
		IMPLICIT NONE
		Integer  :: sunit,i,j,k,l,m,isp, mark, idim
		INTEGER, SAVE :: zone_count = 0
		LOGICAL, SAVE :: first_pass=.TRUE.
		REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
			& mean_vel(ndim), mean_force(ndim), position(ndim)
		CHARaCTER*80 :: FILENAME1, filename2, filename3
		integer, save :: sphrunit1, sphrunit2, sphrunit3
		CHARACTER(LEN=80) :: formfile

		real(prcn), allocatable :: out_arr(:,:,:,:), trans_buf(:)
		integer :: node_num, iproc
		integer :: jj, j1, j2, j3
		real(prcn) :: tmp, tmp1, tmp2, tmp3, mean_energy
		integer :: iphs, part_start, part_end
		logical :: filexist

		real(prcn), allocatable :: acc_fluc(:,:), vel_fluc(:,:), acc_var(:,:), vel_var(:,:), acc_fluc_meanf(:,:), acc_var_meanf(:,:)

		formfile='formatted' 
		j1=1
		j2=my/2
		j3=my

		if (I_AM_NODE_ZERO) write (*,*) "IN FLOW_SNAPSHOT"

#if PARALLEL
		if (I_AM_NODE_ZERO) then
			allocate(out_arr(mx1,3,mz,ndim+3))

			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				! velocity fluctuations for node zero
				do idim=1, ndim+1
					do i=1, nx
						do k=1, mz
							if (idim<=ndim) then
								out_arr(i,jj,k,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
		!							if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
							else
								out_arr(i,jj,k,ndim+1) = pbcp(i,j,k) !-ufmean(idim)
							endif
						enddo
					enddo
				enddo
			enddo

			do iproc=1,nproc-1
				node_num = mz*(ends(iproc)-starts(iproc)+1)*(ndim+1)
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
					call mpi_recv(trans_buf(1), node_num, mpi_double_precision, iproc, iproc, decomp_group, status, err_code)

					l=0
					do idim=1, ndim+1
						do k=1, mz
							do i=starts(iproc),ends(iproc)
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
			enddo
		else
			! recieving velocity fluctuations from node zero
			node_num=mz*nx*(ndim+1)
			allocate(trans_buf(node_num))

			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				l=0
				do idim=1,ndim+1
					do k=1, mz
						do i=1, nx
							l=l+1
							if (idim<=ndim) then
								trans_buf(l) = ubcp(i,j,k,idim) !-ufmean(idim)
	!							if (.not.fluid_atijk(i,j,k)) trans_buf(l) = 0d0
							else
								trans_buf(l) = ubcp(i,j,k,ndim+1) !-ufmean(idim)
							endif
						enddo
					enddo
				enddo

				call mpi_send(trans_buf(1), node_num, mpi_double_precision, node_zero, myid, decomp_group, err_code)
			enddo
			deallocate(trans_buf)
		endif
#else
		allocate(out_arr(mx1,3,mz,ndim+3))
		do idim=1, ndim+1
			do k=1, mz
				do jj=1, 3
					if (jj==1) then
						j=j1
					elseif (jj==2) then
						j = j2
					elseif (jj==3) then
						j = j3
					endif

					do i=1, mx1
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
		write (*,*) "mean_energy = ", mean_energy
		if (I_AM_NODE_ZERO) then
			write (*,*) "GENERATING THE SNAPSHOT OF THE FIELD"

			do k=1, mz
				do j=1, 3
					do i=1, mx1
!						out_arr(i,j,k,1:ndim) = out_arr(i,j,k,1:ndim) / umeanslip
!						out_arr(i,j,k,ndim+1) = out_arr(i,j,k,ndim+1) / mean_energy
						out_arr(i,j,k,ndim+2) = abs(dot_product(out_arr(i,j,k,1:ndim),out_arr(i,j,k,1:ndim)))

						tmp = 0d0
						do idim=1, ndim
							tmp = tmp + (out_arr(i,j,k,idim)-ufmean(idim)) * (out_arr(i,j,k,idim)-ufmean(idim))
						enddo
						out_arr(i,j,k,ndim+3) = tmp / 2 !umeanslip**2
					enddo
				enddo
			enddo

			sunit     = 30
			sphrunit1 = 31 !getnewunit(minunitno,maxunitno)
			sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
			sphrunit3 = 33 !getnewunit(minunitno,maxunitno)

			FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'
			filename2 = TRIM(RUN_NAME)//'_sphr_motion_act1.dat'
			filename3 = TRIM(RUN_NAME)//'_sphr_motion_act2.dat'


			inquire (file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='replace')

				OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='replace')
			ELSE
				OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='old', position='append')

				OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='old', position='append')

			endif

			write(sunit,*)'ZONE T = "', t/t_conv/maxvolfrac, '",'
			write(sunit,*)'DATAPACKING=POINT, I =', mx1,  ', J=', 3, ', K=', mz

			do k=1, mz
				do jj=1, 3
					if (jj==1) then
						j=j1
					elseif (jj==2) then
						j = j2
					elseif (jj==3) then
						j = j3
					endif

					do i=1, mx1
						write(sunit,"(3i6,3d15.7)") i, j, k, out_arr(i,jj,k,4) / mean_energy, out_arr(i,jj,k,5) / umeanslip, out_arr(i,jj,k,6) / mean_energy
					enddo
				enddo
			enddo
			close(sunit,status='keep')

			allocate(acc_fluc(nbody,ndim), vel_fluc(nbody,ndim), acc_var(nphases,ndim), vel_var(nphases,ndim))
			allocate(acc_fluc_meanf(nbody,ndim), acc_var_meanf(nphases,ndim))
			call AiVi(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)

			WRITE(sphrunit1,*)'ZONE T= "', t/t_conv/maxvolfrac, ' " '

			do iphs=1, nphases
				part_start = phase_array(iphs)%pstart
				part_end   = phase_array(iphs)%pend

				do m=part_start, part_end
					tmp1 = sqrt(dot_product(velbdy(m,:),velbdy(m,:))) / umeanslip
					tmp2 = dot_product(velbdy(m,:)-phase_array(iphs)%mean_spec_vel(:), velbdy(m,:)-phase_array(iphs)%mean_spec_vel(:)) / 2 / mean_energy

					if (dot_product(acc_var(iphs,:),vel_var(iphs,:))>post_small) then
						tmp3 = dot_product(acc_fluc(m,:),vel_fluc(m,:)) / dot_product(acc_var(iphs,:),vel_var(iphs,:))
					else
						tmp3 = dot_product(acc_fluc(m,:),vel_fluc(m,:))
					endif

					write (sphrunit1,"(4d15.7)")  xc(m,:), radbdy(m)
!					write (sphrunit1,"(7d15.7)")  xc(m,:), radbdy(m), tmp1, tmp2, tmp3
				enddo
			enddo
			close(sphrunit1)

			deallocate(out_arr)
		endif
	END SUBROUTINE flow_snapshot2



  SUBROUTINE flow_snapshot3
    Use nlarrays , only : ur1, uf1
	 use nlmainarrays, only : ubcp, pbcp
    USE bcsetarrays, ONLY :  omega => fr
    USE dem_mod, only : is_mobile, des_pos_new, des_radius
    IMPLICIT NONE
    Integer  :: sunit,i,j,k,l,m,isp, mark, idim
    INTEGER, SAVE :: zone_count = 0
    LOGICAL, SAVE :: first_pass=.TRUE.
    REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
         & mean_vel(ndim), mean_force(ndim), position(ndim)
    CHARaCTER*80 :: FILENAME1, filename2, filename3
    integer, save :: sphrunit1, sphrunit2, sphrunit3
    CHARACTER(LEN=80) :: formfile

	real(prcn), allocatable :: out_arr(:,:,:,:), trans_buf(:)
	integer :: node_num, iproc
	integer :: kk, k1, k2, k3
	real(prcn) :: tmp, tmp1, tmp2, tmp3, mean_energy
	integer :: iphs, part_start, part_end
	logical :: filexist
	

	formfile='formatted' 

	call compute_omega

	k1=1
	k2=mz/2
	k3=mz

	if (I_AM_NODE_ZERO) write (*,*) "IN FLOW_SNAPSHOT"

#if PARALLEL
	if (I_AM_NODE_ZERO) then
		allocate(out_arr(mx1,my,3,2*ndim+1))

		do kk=1, 3
			if (kk==1) then
				k=k1
			elseif (kk==2) then
				k = k2
			elseif (kk==3) then
				k = k3
			endif

			! velocity fluctuations for node zero
			do idim=1, 2*ndim+1
				do j=1, my
					do i=1, nx

						if (idim<=ndim) then
							out_arr(i,j,kk,idim) = ubcp(i,j,k,idim)
						elseif (idim<=2*ndim) then
							out_arr(i,j,kk,idim) = omega(i,j,k,idim-ndim)
						else
							out_arr(i,j,kk,idim) = pbcp(i,j,k)
						endif

					enddo
				enddo
			enddo
		enddo

		do iproc=1,nproc-1
			node_num = my*(ends(iproc)-starts(iproc)+1)*(2*ndim+1)
			allocate(trans_buf(node_num))

			do kk=1, 3
				if (kk==1) then
					k = k1
				elseif (kk==2) then
					k = k2
				elseif (kk==3) then
					k = k3
				endif

				! collecting velocity fluctuations from other processes
				call mpi_recv(trans_buf(1), node_num, mpi_double_precision, iproc, iproc, decomp_group, status, err_code)

				l=0
				do idim=1, 2*ndim+1
					do j=1, my
						do i=starts(iproc),ends(iproc)
							l=l+1

							out_arr(i,j,kk,idim) = trans_buf(l)
						enddo
					enddo
				enddo
			enddo
			deallocate(trans_buf)
		enddo
	else
		! recieving velocity fluctuations from node zero
		node_num=my*nx*(2*ndim+1)
		allocate(trans_buf(node_num))

		do kk=1, 3
			if (kk==1) then
				k=k1
			elseif (kk==2) then
				k = k2
			elseif (kk==3) then
				k = k3
			endif

			l=0
			do idim=1, 2*ndim+1
				do j=1, my
					do i=1, nx
						l=l+1
						if (idim<=ndim) then
							trans_buf(l) = ubcp(i,j,k,idim)
						elseif (idim<=2*ndim) then
							trans_buf(l) = omega(i,j,k,idim-ndim)
						else
							trans_buf(l) = pbcp(i,j,k)
						endif
					enddo
				enddo
			enddo

			call mpi_send(trans_buf(1),node_num,mpi_double_precision,node_zero,myid,decomp_group,err_code)
		enddo
		deallocate(trans_buf)
	endif
#else
	allocate(out_arr(mx1,my,3,2*ndim+1))

	do idim=1, 2*ndim+1
		do j=1, my
			do kk=1, 3
				if (kk==1) then
					k=k1
				elseif (kk==2) then
					k = k2
				elseif (kk==3) then
					k = k3
				endif

				do i=1, mx1
					if (idim<=ndim) then
						out_arr(i,j,kk,idim) = ubcp(i,j,k,idim)
					elseif (idim<=2*ndim) then
						out_arr(i,j,kk,idim) = omega(i,j,k,idim-ndim)
					else
						out_arr(i,j,kk,idim) = pbcp(i,j,k)
					endif
				enddo
			enddo
		enddo
	enddo
#endif

	mean_energy = 0.5*umeanslip**2
	if (I_AM_NODE_ZERO) then
		write (*,*) "GENERATING THE SNAPSHOT OF THE FIELD"

		out_arr(:,:,:,1:ndim) = out_arr(:,:,:,1:ndim) / umeanslip

		out_arr(:,:,:,ndim+1:2*ndim) = out_arr(:,:,:,ndim+1:2*ndim) / (umeanslip/dia_phys)
		out_arr(:,:,:,2*ndim+1) = out_arr(:,:,:,2*ndim+1) / mean_energy

      sunit     = 30
		sphrunit1 = 31 !getnewunit(minunitno,maxunitno)
		sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
		sphrunit3 = 33 !getnewunit(minunitno,maxunitno)

		FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'

		inquire (file=trim(filename1), exist=filexist)
		if (.not.filexist) then
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='replace')
			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='unknown')
		ELSE
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='old', position='append')
			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='old', position='append')
		endif

		write(sunit,*)'ZONE T = "', t/t_conv/(1-maxvolfrac), '",'
		write(sunit,*)'DATAPACKING=POINT, I =', mx1,  ', J=', my, ', K=', 3

		do kk=1, 3
			if (kk==1) then
				k=k1
			elseif (kk==2) then
				k = k2
			elseif (kk==3) then
				k = k3
			endif

			do j=1, my
				do i=1, mx1
					write(sunit,"(3i6,7d15.7)") i, j, k, out_arr(i,j,kk,:) !out_arr(i,jj,k,1), out_arr(i,jj,k,2), out_arr(i,jj,k,3), out_arr(i,jj,k,4)
				enddo
			enddo
		enddo
		close(sunit,status='keep')

		WRITE(sphrunit1,*)'ZONE T= "', t/t_conv/(1-maxvolfrac), ' " '

		do iphs=1, nphases
			part_start = phase_array(iphs)%pstart
			part_end   = phase_array(iphs)%pend

			do m=part_start, part_end
				write (sphrunit1,"(5d15.7)")  xc(m,:), radbdy(m), sqrt(dot_product(velbdy(m,:),velbdy(m,:)))/umeanslip
			enddo
		enddo
		close(sphrunit1)
		deallocate(out_arr)
	endif
  END SUBROUTINE flow_snapshot3



	subroutine write_drag
		implicit none
		character*50 filename1
		integer :: line, nvar, nvar1, nvar2, iphs
		real(prcn) :: tmp1, tmp2, tmp3, tmp4, ibm_drag

		if (I_AM_NODE_ZERO) then
			if (.not. post_no_flow_mem_alloc) then
				filename1=trim(run_name)//"_drag.dat"
				open (unit=1, file=trim(filename1), status="replace", action="write")
				write (1,*) "zone"

				if (nphases==1) then
					call compute_ibm_drag(maxvolfrac, re, ibm_drag)

					write (1,"(4f10.4,4D15.7)") dbydx, lybyd, re, maxvolfrac, &
											& 	norm_drag, norm_drag_spec(1), ibm_drag, abs(norm_drag_spec(1)-ibm_drag) / ibm_drag
				else
					write (1,"(8f10.4,3D15.7)") dbydx, lybyd, re, maxvolfrac, &
											& 	(phase_array(iphs)%volfracg, iphs=1,nphases), &
											&	(yalpha(iphs),iphs=1,nphases), 			&
											&	norm_drag, norm_drag_spec(1:nphases)
				endif
				close (1)

				filename1=trim(run_name)//"_drag_parts.dat"
				open (unit=1, file=trim(filename1), status="replace", action="write")
				write (1,*) "zone"
				tmp1 = sqrt(dot_product(pres_total,pres_total)) / nbody
				tmp2 = sqrt(dot_product(visc_total,visc_total)) / nbody
				tmp3 = sqrt(dot_product(mpg,mpg)) * pi*dia_phys**3/6
				tmp4 = 3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*dia_phys

				write (1,"(4f10.4,3D15.7)") dbydx, lybyd, re, maxvolfrac, &
											& 	(tmp1+tmp3)/tmp4, tmp2/tmp4, (tmp1+tmp2+tmp3)/tmp4
				close(1)
			else
				if (nphases==1) then
					nvar  = 8
					nvar1 = 4
					nvar2 = 4
				else
					nvar  = 11
					nvar1 = 8
					nvar2 = 3
				endif

				filename1 = "_drag.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line)

				nvar  = 7
				nvar1 = 4
				nvar2 = 3

				filename1 = "_drag_parts.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, line)
			endif
		endif
	end subroutine write_drag

	subroutine compute_interphase_transfer
		implicit none

		integer :: istep, jstep, kstep, nstep1, nstep2, io, nvar1, nvar2, nvar, nline
		integer :: max_step = 1000000

		logical :: filexist, match
		character*50 filename1

		real(prcn) :: tmp1, tmp2, x1, x2
		real(prcn), allocatable :: var1(:,:), var2(:,:), var1_tmp(:), var2_tmp(:)
		real(prcn), allocatable :: var_u(:,:), var_d(:,:), var_out(:,:), var_tmp(:)
		integer :: des_var1, des_var2, nstep_tmp, count


		call screen_separator(30,'^')
		write (*,*) "IN COMPUTE_INTERPHASE_TRANSFER"		

		if (.not.post_no_flow_mem_alloc) then
			if (nphases<2) then
				write (*,*) " NOT FOR MONODISPERSE SUSPENSIONS, EXIT..."
				return
			endif

			call compute_mix_mean_vel

			!^^^ reading vel_info.dat ^^^^
			filename1=trim(run_name)//"_vel_info.dat"
			inquire(file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				write (*,*) 'FILE "'//trim(filename1)//' DOES NOT EXIST.'
				stop
			endif
			open (unit=1, file=trim(filename1), status="old", action="read")

			if (nphases==1) then
				nvar1 = 12
			elseif (nphases==2) then
				nvar1 = 16
			endif

			allocate(var1(nvar1,max_step), var1_tmp(nvar1))
			var1_tmp = zero
			var1 = zero

			des_var1 = 1
			nstep1 = 0
			count = 0
			do
				count = count+1
				if (nphases==1) then
					read(1,"(12(2x,g17.8))",iostat=io) var1_tmp (1:nvar1)
				elseif (nphases==2) then
					read(1,"(16(2x,g17.8))",iostat=io) var1_tmp (1:nvar1)
				endif

				if (io>0) then
					write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN LINE ", nstep1+1
				elseif (io==0) then
					if (nstep1>1 .and. var1_tmp(des_var1)<=var1(des_var1,nstep1)) then

!						write (*,*) "a problem in vel : ", count
!						read (*,*)

						do nstep_tmp=nstep1, 1, -1
							if (var1_tmp(des_var1)>var1(des_var1,nstep_tmp)) exit
						enddo
						nstep_tmp = nstep_tmp+1

						var1(:,nstep_tmp:nstep1) = zero
						nstep1 = nstep_tmp
					else
						nstep1 = nstep1 + 1
					endif

					var1(1:nvar1,nstep1) = var1_tmp (1:nvar1)
				elseif (io<0) then
					write (*,*) "END OF FILE ", trim(filename1)
					write (*,*) "NUMBER OF INPUTS = ", nstep1
					exit
				endif
			enddo
			close (1)
			!---------------------------

			!^^^ reading norm_drag.dat ^^^^
			filename1=trim(run_name)//"_norm_drag.dat"
			inquire(file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				write (*,*) 'FILE "'//trim(filename1)//' DOES NOT EXIST.'
				stop
			endif
			open (unit=1, file=trim(filename1), status="old", action="read")

			if (nphases==1) then
				nvar2 = 11
			elseif (nphases==2) then
				nvar2 = 12
			endif

			allocate(var2(nvar2,max_step), var2_tmp(nvar2))
			var2_tmp = zero
			var2 = zero

			des_var2 = 1
			nstep2 = 0
			count = 0
			do
				count = count+1
				if (nphases==1) then
					read(1,"(11(2x, e20.12))",iostat=io) var2_tmp (1:nvar2)
				elseif (nphases==2) then
					read(1,"(12(2x, e20.12))",iostat=io) var2_tmp (1:nvar2)
				endif

				if (io>0) then
					write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN LINE ", nstep2+1
				elseif (io==0) then
					if (nstep2>1 .and. var2_tmp(des_var2)<=var2(des_var2,nstep1)) then

!						write (*,*) "a problem in drag : ", count
!						read (*,*)

						do nstep_tmp=nstep2, 1, -1
							if (var2_tmp(des_var2)>var2(des_var2,nstep_tmp)) exit
						enddo
						nstep_tmp = nstep_tmp+1

						var2(:,nstep_tmp:nstep2) = zero
						nstep2 = nstep_tmp
					else
						nstep2 = nstep2 + 1
					endif

					var2(1:nvar2,nstep2) = var2_tmp (1:nvar2)
				elseif (io<0) then
					write (*,*) "END OF FILE ", trim(filename1)
					write (*,*) "NUMBER OF INPUTS = ", nstep2
					exit
				endif
			enddo
			close (1)
			!---------------------------

			allocate(var_u(10,nstep1))
			allocate(var_d(6,nstep2))
			allocate(var_out(29,nstep1))

			do istep=1, nstep1
				x1 = phase_array(1)%volfrac !/maxvolfrac ! volfrac_rat
				x2 = phase_array(2)%volfrac !/maxvolfrac ! volfrac_rat

				var_u(1,istep) = var1(2,istep) / (one-maxvolfrac) !time

				mix_mean_slip_mag = var1(4,istep) * vis / (1-maxvolfrac) / char_length

				mix_mean_vel(:) = var1(5:7,istep) * mix_mean_slip_mag

				mix_mean_vel_mag = sqrt(dot_product(mix_mean_vel,mix_mean_vel))

				ufmean_mag = mix_mean_slip_mag+mix_mean_vel_mag

				var_u(2,istep) = ufmean_mag
				var_u(3,istep) = mix_mean_vel_mag
				var_u(4,istep) = mix_mean_slip_mag

				phase_array(1)%mean_spec_vel(1:ndim) = (var1(5:7,istep)-x2*var1(8:10,istep))/(x1+x2) 
				phase_array(2)%mean_spec_vel(1:ndim) = (var1(5:7,istep)+x1*var1(8:10,istep))/(x1+x2)

				var_u(5,istep) = sqrt(dot_product(phase_array(1)%mean_spec_vel,phase_array(1)%mean_spec_vel)) * mix_mean_slip_mag
				var_u(6,istep) = sqrt(dot_product(phase_array(2)%mean_spec_vel,phase_array(2)%mean_spec_vel)) * mix_mean_slip_mag

				var_u(7,istep) = (var1(11,istep)*mix_mean_slip_mag)**2
				var_u(8,istep) = (var1(12,istep)*mix_mean_slip_mag)**2

				var_u(9, istep) = rhos * phase_array(1)%volfrac * (var1(11,istep)*mix_mean_slip_mag)**2 * three/two
				var_u(10,istep) = rhos * phase_array(2)%volfrac * (var1(12,istep)*mix_mean_slip_mag)**2 * three/two
			enddo

			do istep=1, nstep2
				var_d(1,istep) = var2(1,istep) / (one-maxvolfrac)

				var_d(2,istep) = var2(4,istep)

				var_d(3,istep) = var2(5,istep)
				var_d(4,istep) = var2(6,istep)

				var_d(5,istep) = var2(5,istep) *3*pi*vis*phase_array(1)%dia*(one-maxvolfrac)*mix_mean_slip_mag * phase_array(1)%npart
				var_d(6,istep) = var2(6,istep) *3*pi*vis*phase_array(2)%dia*(one-maxvolfrac)*mix_mean_slip_mag * phase_array(2)%npart
			enddo

			istep = 0
			jstep = 0
			kstep = 0
			do istep=1, nstep1
				do
					match = .false.
					jstep = jstep+1
					if (jstep>nstep2) then
						jstep = 0
						exit
					endif

					if (abs(var_u(1,istep)-var_d(1,jstep))<post_small) then
						match = .true.
						exit
					endif
				enddo

				if (match==.true.) then
					kstep = kstep + 1
					var_out(1,kstep) = var_u(1,istep) !TIME

					ufmean_mag = var_u(2,istep)
					mix_mean_vel_mag = var_u(3,istep)
					mix_mean_slip_mag = var_u(4,istep)

					usmean1_mag = var_u(5,istep)
					usmean2_mag = var_u(6,istep)

					var_out(2,kstep) = var_u(2,istep)/mix_mean_slip_mag !UFMEAN/MEANSLIP
					var_out(3,kstep) = var_u(3,istep)/mix_mean_slip_mag !USMEAN/MEANSLIP

					var_out(4,kstep) = usmean1_mag/mix_mean_slip_mag
					var_out(5,kstep) = usmean2_mag/mix_mean_slip_mag

					var_out(6,kstep) = sqrt(dot_product(var1(8:10,istep),var1(8:10,istep))) * mix_mean_slip_mag
					var_out(7,kstep) = usmean2_mag

					! ^^^^ MEAN DRAG ^^^^
					var_out(8,kstep) = var_d(3,jstep) !NORMALIZED FORCE PER PARTICLE
					var_out(9,kstep) = var_d(4,jstep) !NORMALIZED FORCE PER PARTICLE

					var_out(10,kstep) = var_d(3,jstep) * phase_array(1)%npart / (doml(1)*doml(2)*doml(3)) !NORMALIZED FORCE PER VOLUME
					var_out(11,kstep) = var_d(4,jstep) * phase_array(1)%npart / (doml(1)*doml(2)*doml(3)) !NORMALIZED FORCE PER VOLUME

					var_out(12,kstep) = var_d(5,jstep) / (doml(1)*doml(2)*doml(3)) !TOTAL FORCE PER VOLUME / IF1
					var_out(13,kstep) = var_d(6,jstep) / (doml(1)*doml(2)*doml(3)) !TOTAL FORCE PER VOLUME / IF2

					! ^^^^ MEAN TKE TRANSFER ^^^^
					var_out(14,kstep) = abs(usmean1_mag-ufmean_mag) * var_d(4,jstep) / (vis*(mix_mean_slip_mag/dia_phys)**2) &
									& / ((one-maxvolfrac)*doml(1)*doml(2)*doml(3))
					var_out(15,kstep) = abs(usmean2_mag-ufmean_mag) * var_d(5,jstep) / (vis*(mix_mean_slip_mag/dia_phys)**2) &
									& / ((one-maxvolfrac)*doml(1)*doml(2)*doml(3))

					var_out(16,kstep) = var_out(12,kstep)+var_out(13,kstep)
				
					var_out(17,kstep) = abs(usmean1_mag-ufmean_mag) * var_d(5,jstep) / (vis*(mix_mean_slip_mag/dia_phys)**2) &
									& / phase_array(1)%npart
					var_out(18,kstep) = abs(usmean2_mag-ufmean_mag) * var_d(6,jstep) / (vis*(mix_mean_slip_mag/dia_phys)**2) &
									& / phase_array(2)%npart

					var_out(19,kstep) = var_u(7,istep)/var_u(8,istep)
					var_out(20,kstep) = var_u(9,istep)/var_u(10,istep)
				else
					write (*,*) "ISTEP DIDNOT MATCH : ", istep, kstep
					read (*,*)
				endif

			enddo

			do istep=1, kstep
				if (istep>1) then
					dt = var_out(1,istep)-var_out(1,istep-1)
					var_out(21,istep) = (var_out(6,istep)-var_out(6,istep-1))/dt
				else
					var_out(21,istep) = 0
				endif
			enddo

			allocate(var_tmp(kstep))
			var_tmp = zero
			do istep=1, kstep
				if (istep<=skip_num) then
					var_tmp(istep) = sum(var_out(21,1:2*istep-1)) / (2*istep-1)
				elseif (istep>kstep-skip_num) then
					var_tmp(istep) = sum(var_out(21,2*istep-kstep:kstep)) / (2*kstep-2*istep+1)
				else
					var_tmp(istep) = sum(var_out(21,istep-skip_num:istep+skip_num)) / (2*skip_num+1)
				endif
			enddo
			var_out(21,:) = var_tmp(:)

			var_tmp = zero
			do istep=1, kstep
				if (istep<=skip_num) then
					var_tmp(istep) = sum(var_out(12,1:2*istep-1)) / (2*istep-1)
				elseif (istep>kstep-skip_num) then
					var_tmp(istep) = sum(var_out(12,2*istep-kstep:kstep)) / (2*kstep-2*istep+1)
				else
					var_tmp(istep) = sum(var_out(12,istep-skip_num:istep+skip_num)) / (2*skip_num+1)
				endif
			enddo
			var_out(22,:) = var_tmp(:)

			var_tmp = zero
			do istep=1, kstep
				if (istep<=skip_num) then
					var_tmp(istep) = sum(var_out(13,1:2*istep-1)) / (2*istep-1)
				elseif (istep>kstep-skip_num) then
					var_tmp(istep) = sum(var_out(13,2*istep-kstep:kstep)) / (2*kstep-2*istep+1)
				else
					var_tmp(istep) = sum(var_out(13,istep-skip_num:istep+skip_num)) / (2*skip_num+1)
				endif
			enddo
			var_out(23,:) = var_tmp(:)

			do istep=1, kstep
				var_out(24,istep) = var_out(12,istep) / rhos/phase_array(1)%volfrac
				var_out(25,istep) = var_out(13,istep) / rhos/phase_array(2)%volfrac

				var_out(26,istep) = var_out(22,istep) / rhos/phase_array(1)%volfrac
				var_out(27,istep) = var_out(23,istep) / rhos/phase_array(2)%volfrac

				var_out(28,istep) = -((var_out(21,istep) - (var_out(24,istep) - var_out(25,istep)))) / (1/rhos/phase_array(1)%volfrac+1/rhos/phase_array(2)%volfrac) * dia_phys**2/18/vis/maxvolfrac/(1-maxvolfrac)/mix_mean_slip_mag 
				var_out(29,istep) = -((var_out(21,istep) - (var_out(26,istep) - var_out(27,istep)))) / (1/rhos/phase_array(1)%volfrac+1/rhos/phase_array(2)%volfrac) * dia_phys**2/18/vis/maxvolfrac/(1-maxvolfrac)/mix_mean_slip_mag 

			enddo

			filename1=trim(run_name)//"_interphase_transfer.dat"
			open (unit=1, file=trim(filename1), status="replace", action="write")
			write (1,*) "zone"

			do istep=2, kstep
				write (1,"(1f8.4,28d15.7)") var_out(:,istep)
			enddo
			close(1)

			deallocate(var_u, var_d, var_out, var1, var2, var1_tmp, var2_tmp)
		else
			if (post_no_flow_mem_alloc) then
				nvar  = 29
				nvar1 = 1
				nvar2 = 28

				line = nbins

				filename1 = "_interphase_transfer.dat"
				call mis_average(nvar, nvar1, nvar2, filename1, nbins)
			endif
		endif

		call screen_separator(30,'^')
	end subroutine compute_interphase_transfer


	subroutine fluid_particle_acceleration
!		use boundary_condition
		use bcsetarrays, only : ppr, diffn
		implicit none

		real(prcn), allocatable :: f_acc(:,:), p_acc(:,:), wt(:), hist(:), radbin(:), tmp_array(:)
		real(prcn) :: u_min, u_max
		integer :: ibody, idim, i, j, k, count, unit1, nvar, nvar1, nvar2, ip, im
		character*50 filename1, filename2

		call screen_separator(30,'^')
		write (*,*) "IN FLUID_PARTICLE_ACCELERATION"

		if (.not.post_no_flow_mem_alloc) then
			allocate(f_acc(ndim, count_fluid), p_acc(ndim, count_fluid))

			f_acc = zero
			p_acc = zero

!			call calc_pgrad
!			call calc_visc


			count = 0
			do k=1, mz
				do j=1, my
					do i=1, nx
						if (fluid_atijk(i,j,k)) then
							count = count+1
							im = i-1
							if (im<1) im = nx

							ip = i

							do idim=1, ndim
								f_acc(idim, count) = diffn(i,j,k, idim) + half * (ppr(im,j,k,idim)+ppr(ip,j,k,idim))
							enddo
						endif
					enddo
				enddo
			enddo

			do ibody=1, nbody
				do idim=1, ndim
					p_acc(idim, ibody) = visc(idim, ibody) + pres(idim, ibody)
				enddo
			enddo

			unit1 = 1
			filename1 = trim(run_name)//"_pdf_f_acc.dat"
			open  (unit=unit1, file=trim(filename1), status="replace")

			allocate(wt(count_fluid))
			wt = 1d0/count_fluid


			allocate(hist(nbins), radbin(nbins))
			allocate(tmp_array(count_fluid))


			do idim=1, ndim
!				tmp_array(:) = f_acc(idim,:)

if (idim==1) then
	filename2 = trim(run_name)//"_ravi.dat"
	j=0
	open (2, file=trim(filename2), status="replace")
	do i=1, count_fluid
		if (mod(i,200)==0) then
			j=j+1
			tmp_array(j) = f_acc(idim,i)
			write (2,"(1e15.7)") tmp_array(j)
		endif
	enddo
	close (2)
endif



				call make_histogram(tmp_array(1:j), j, nbins, radbin, hist)

				write (unit1,*) "zone"
				do i=1, nbins
					if (hist(i)>small_number) write (1,"(2d15.7)") radbin(i), hist(i)
				enddo
			enddo
			close (unit1)

			unit1 = 1
			filename1 = trim(run_name)//"_pdf_p_acc.dat"
			open  (unit=unit1, file=trim(filename1), status="replace")

			deallocate(tmp_array)
			allocate(tmp_array(nbody))
			do idim=1, ndim
				tmp_array(:) = p_acc(idim,:)
				call make_histogram(tmp_array, nbody, nbins, radbin, hist)

				write (unit1,*) "zone"
				do i=1, nbins
					if (hist(i)>small_number) write (1,"(2d15.7)") radbin(i), hist(i)
				enddo
			enddo
			close (unit1)

		else
			nvar  = 2
			nvar1 = 1
			nvar2 = 1

			line = nbins

			filename1 = "_pdf_f_acc.dat"
			call mis_average(nvar, nvar1, nvar2, filename1, nbins)

			filename1 = "_pdf_p_acc.dat"
			call mis_average(nvar, nvar1, nvar2, filename1, nbins)
		endif
		call screen_separator(30,'-')

	end subroutine fluid_particle_acceleration


	subroutine make_histogram(array, n, nbins, radbin, hist)
		implicit none
		integer, intent(in) :: n, nbins
		real(prcn), intent(inout) :: array(n)
		real(prcn), intent(out) :: radbin(nbins), hist(nbins)

!		real(prcn) :: tmp_array(n)
		real(prcn) :: mean, var, sd, left, right, dr
		integer :: i, ibin

		call calc_avr_var(n, array, mean, var)

		sd = sqrt(var/n)

		array(:) = (array(:) - mean) / sd

write (*,*) mean, sd

		call calc_avr_var(n, array, mean, var)

		sd = sqrt(var/n)

write (*,*) mean, sd

		left = minval(array(:))
		right = maxval(array(:))

		write (*,*) "LEFT, RIGHT = ", left, right

		dr = (right-left) / nbins

		do i=1, nbins
			radbin(i) = left + (i-.5)*dr
		enddo

		hist = zero
		do i=1, n
			ibin = (array(i)-left) / dr + 1
			if (ibin>nbins) ibin = nbins
			hist(ibin) = hist(ibin) + 1
		enddo

		hist(:) = hist(:) / n / dr

		write (*,*) "SUM OF HIST = ", sum(hist(:)) * dr

		radbin(:) = radbin(:)
	end subroutine make_histogram


	subroutine compute_gofr_avg
		implicit none
		character*50 filename
		integer :: i
		integer :: nvar, nvar1, nvar2

		if (post_no_flow_mem_alloc) then
			nvar  = 3
			nvar1 = 1
			nvar2 = 2

			line = nbins

			filename = "_gofr_3d.dat"
			
			call mis_average(nvar, nvar1, nvar2, filename, line)
		endif
	end subroutine compute_gofr_avg

	subroutine velocity_output
		implicit none
		character*50 filename1
		integer :: i,j,k, idim
		real(prcn) :: tmp, mean_energy
		
		filename1 = trim(run_name)//"_velcity_field.dat"
		open (unit=1, file=trim(filename1), status="replace")
!		write (1,*) "variables = x,y,z,u,v,w"
		write (1,*) "zone"
		write (1,"(3(1A,1I),1a)") "i=", nx, " j=", my, " k=", mz, " f=point"


		call compute_mix_mean_vel

		mean_energy = half*mix_mean_slip_mag**2

		do k=1, mz
			do j=1, my
!				do i=nx/2+1, nx
!					if (fluid_atijk(i,j,k)) then
!						write (1,"(3i6,3d15.7)") i-(nx+1), j, k, ubcp(i,j,k,:)
!					else
!						write (1,"(3i6,3d15.7)") i-(nx+1), j, k, zero, zero, zero
!					endif
!				enddo

				do i=1,nx
!					if (fluid_atijk(i,j,k)) then
						tmp = zero
						do idim=1,ndim
							tmp = tmp + (ubcp(i,j,k,idim)-ufmean(idim)) * (ubcp(i,j,k,idim)-ufmean(idim)) / 2
						enddo
						write (1,"(3i6,1d15.7)") i, j, k, tmp / mean_energy
!					else
!						write (1,"(3i6,3d15.7)") i, j, k, zero, zero, zero
!					endif
				enddo
			enddo
		enddo

		close (1)
	end subroutine velocity_output

	subroutine volumetric_drag
		use bcsetarrays
		use nlmainarrays
		implicit none
		integer    :: idelta, i, j, k, ib, ie, jb, je, kb, ke, ii, jj, kk, iii, jjj, kkk, dim1, dim2, m, ibody, idim
		real(prcn) :: dx_filter, dy_filter, dz_filter, rmax, vol_drag_mag, norm_force, rad

		real(prcn), allocatable :: vol_drag(:,:)
		real(prcn) :: box_center(ndim), box_limit(ndim,2), dist(ndim), vol_drag_average(ndim)

		character*50 filename1

		integer :: nx_new, my_new, mz_new, vol_count, vol_num
		logical :: inbox

		call screen_separator(30,'^')
		write (*,*) "IN VOLUMETRIC_DRAG..."

		if (.not.post_no_flow_mem_alloc) then


			call compute_mix_mean_vel

			norm_force = 3*pi*dia_phys*vis*mix_mean_slip_mag*(1-maxvolfrac)

			filename1 = trim(run_name)//"_volumetric_drag.dat"
			open (unit=1, file=trim(filename1), status="replace")
			write (1,*) "zone"


			vol_num = nx*my*mz
			allocate(vol_drag(vol_num,ndim))


write (*,*) "NBODY = ", nbody
write (*,*) "NBND = ", nbnd
write (*,*) "VOL_COUNT = ", vol_count

			do idelta=1, nx, 2
				write (*,*) "DELTA = ", idelta

				nx_new = nx / idelta
				my_new = my / idelta
				mz_new = mz / idelta

				dx_filter = idelta*dx
				dy_filter = idelta*dy
				dz_filter = idelta*dz

				vol_drag = zero

				vol_count = 0
				do k=1, mz
					do j=1, my
						do i=1, nx
							vol_count = vol_count + 1

if (mod(vol_count,nx*my)==1) write (*,*) "VOL_COUNT = ", vol_count

#if 0
							ib = i-idelta/2
							ie = i+idelta/2

							jb = j-idelta/2
							je = j+idelta/2

							kb = k-idelta/2
							ke = k+idelta/2

							do kk=kb, ke
								if (kk<1) then
									kkk = kk + mz
								elseif (kk>mz) then
									kkk = kk - mz
								else
									kkk = kk
								endif

								do jj=jb, je
									if (jj<1) then
										jjj = jj + my
									elseif (jj>my) then
										jjj = jj - my
									else
										jjj = jj
									endif

									do ii=ib, ie
										if (ii<1) then
											iii = ii + nx
										elseif (ii>nx) then
											iii = ii - nx
										else
											iii = ii
										endif

										total_num = total_num + 1
										if (.not.fluid_atijk(iii,jjj,kkk)) solid_num = solid_num+1
									enddo
								enddo
							enddo
#endif					

							box_center(1) = i
							box_center(2) = j
							box_center(3) = k

							box_limit(1,1) = -dble(idelta)/2
							box_limit(1,2) =  dble(idelta)/2
							box_limit(2,1) = -dble(idelta)/2
							box_limit(2,2) =  dble(idelta)/2
							box_limit(3,1) = -dble(idelta)/2
							box_limit(3,2) =  dble(idelta)/2

							rmax = dble(my)/2
							do ibody=1,nbody
								dist(:) = abs(xc(ibody,:) - box_center(:))
								do idim=1, ndim
									if (dist(idim)>rmax) dist(idim) = 2*rmax - dist(idim)
								enddo
								rad = sqrt(dot_product(dist(:),dist(:)))
								if (rad>sqrt(3.)*idelta/2.+radbdy(ibody)) goto 10
!write (*,*) ibody


								do m=1, nbnd
									dist(:) = xc_bnd(ibody,m,:) - box_center(:)

									do idim=1, ndim
										if (dist(idim)>rmax) then
											dist(idim) = dist(idim) - 2*rmax
										elseif (dist(idim)<-rmax) then
											dist(idim) = dist(idim) + 2*rmax
										endif
									enddo

									do idim=1, ndim
										if (box_limit(idim,1)<=dist(idim) .and. dist(idim)<=box_limit(idim,2)) then
											inbox = .true.
										else
											inbox = .false.
											exit
										endif
									enddo

									if (inbox) then
										vol_drag(vol_count,:) = vol_drag(vol_count,:) + (visc_force_bnd(ibody,m,:)+pres_force_bnd(ibody,m,:))
									endif
								enddo
10								continue
							enddo

						enddo
					enddo
				enddo

write (*,*) 1
				if (vol_count/=vol_num) then
					write (*,*) "VOL_NUM/=VOL_COUNT"
					stop
				endif

write (*,*) 2

				vol_drag_average = zero
				do m=1, vol_num
					vol_drag_average(:) = vol_drag_average(:) + vol_drag(m,:)
				enddo

write (*,*) 3

				vol_drag_average(:) = vol_drag_average(:) / vol_num
				vol_drag_average(:) = vol_drag_average(:) - mpg(:) * maxvolfrac * (dble(idelta)/dbydx*dia_phys)**3

				vol_drag_mag = sqrt(dot_product(vol_drag_average(:),vol_drag_average(:)))

write (*,*) 4

				write (1,"(6d15.7)") dble(idelta)/dbydx, (dble(idelta)/my)**3, vol_drag_average(:)/norm_force, vol_drag_mag/norm_force

write (*,*) 5

write (*,*) dble(idelta)/dbydx, vol_drag_mag/norm_force
read (*,*)
			enddo
			close (1)
		else
		endif
	end subroutine volumetric_drag








	subroutine volumetric_drag2
		implicit none
		integer    :: idelta, i, j, k, ibody, idim
		real(prcn) :: rmax, vol_drag_mag, norm_force, rad

		real(prcn), allocatable :: vol_drag(:,:)
		real(prcn) :: box_center(ndim), box_limit(ndim,2), dist(ndim), vol_drag_average(ndim)

		character*50 filename1

		integer :: vol_count, vol_num
		logical :: inbox

		call screen_separator(30,'^')
		write (*,*) "IN VOLUMETRIC_DRAG..."

		if (.not.post_no_flow_mem_alloc) then


			call compute_mix_mean_vel

			norm_force = 3*pi*dia_phys*vis*mix_mean_slip_mag*(1-maxvolfrac)

			filename1 = trim(run_name)//"_volumetric_drag.dat"
			open (unit=1, file=trim(filename1), status="replace")
			write (1,*) "zone"


			vol_num = nx*my*mz
			allocate(vol_drag(vol_num,ndim))


			write (*,*) "NBODY = ", nbody
			!write (*,*) "NBND = ", nbnd
			write (*,*) "VOL_COUNT = ", vol_num

			do idelta=1, nx, 2
				write (*,*) "DELTA = ", idelta

!				nx_new = nx / idelta
!				my_new = my / idelta
!				mz_new = mz / idelta

!				dx_filter = idelta*dx
!				dy_filter = idelta*dy
!				dz_filter = idelta*dz

				vol_count = 0
				do k=1, mz
					do j=1, my
						do i=1, nx
							vol_count = vol_count + 1

if (mod(vol_count,nx*my)==1) write (*,*) "VOL_COUNT = ", vol_count

#if 0
							ib = i-idelta/2
							ie = i+idelta/2

							jb = j-idelta/2
							je = j+idelta/2

							kb = k-idelta/2
							ke = k+idelta/2

							do kk=kb, ke
								if (kk<1) then
									kkk = kk + mz
								elseif (kk>mz) then
									kkk = kk - mz
								else
									kkk = kk
								endif

								do jj=jb, je
									if (jj<1) then
										jjj = jj + my
									elseif (jj>my) then
										jjj = jj - my
									else
										jjj = jj
									endif

									do ii=ib, ie
										if (ii<1) then
											iii = ii + nx
										elseif (ii>nx) then
											iii = ii - nx
										else
											iii = ii
										endif

										total_num = total_num + 1
										if (.not.fluid_atijk(iii,jjj,kkk)) solid_num = solid_num+1
									enddo
								enddo
							enddo
#endif					

							box_center(1) = i
							box_center(2) = j
							box_center(3) = k

							box_limit(1,1) = -dble(idelta)/2
							box_limit(1,2) =  dble(idelta)/2
							box_limit(2,1) = -dble(idelta)/2
							box_limit(2,2) =  dble(idelta)/2
							box_limit(3,1) = -dble(idelta)/2
							box_limit(3,2) =  dble(idelta)/2

							rmax = dble(my)/2
							do ibody=1,nbody
								dist(:) = xc(ibody,:) - box_center(:)
								do idim=1, ndim
									if (dist(idim)>rmax) then
										dist(idim) = dist(idim) - 2*rmax
									elseif (dist(idim)<-rmax) then
										dist(idim) = dist(idim) + 2*rmax
									endif
								enddo

								do idim=1, ndim
									if (box_limit(idim,1)<=dist(idim) .and. dist(idim)<=box_limit(idim,2)) then
										inbox = .true.
									else
										inbox = .false.
										exit
									endif
								enddo

								if (inbox) vol_drag(vol_count,:) = vol_drag(vol_count,:) + force(ibody,:)
							enddo

						enddo
					enddo
				enddo

				if (vol_count/=vol_num) then
					write (*,*) "VOL_NUM/=VOL_COUNT"
					stop
				endif

				vol_drag_average = zero
				do i=1, vol_num
					vol_drag_average(:) = vol_drag_average(:) + vol_drag(i,:)
				enddo

				vol_drag_average(:) = vol_drag_average(:) / vol_num
!				vol_drag_average(:) = vol_drag_average(:) - mpg(:) * maxvolfrac * (dble(idelta)/dbydx*dia_phys)**3

				vol_drag_mag = sqrt(dot_product(vol_drag_average(:),vol_drag_average(:)))

				write (1,"(7d15.7)") dble(idelta)/dbydx, (dble(idelta)/my)**3, vol_drag_average(:)/norm_force, vol_drag_mag/norm_force, vol_drag_mag/norm_force/(dble(idelta)/my)**3
			enddo
			close (1)
		else
		endif
	end subroutine volumetric_drag2




	subroutine filter
		use bcsetarrays
		use nlmainarrays
		implicit none
		integer    :: idelta, i, j, k, ib, ie, jb, je, kb, ke, ii, jj, kk, iii, jjj, kkk, dim1, dim2
		real(prcn) :: tmp_scal1, tmp_scal2, tmp_vec(ndim), meanslip_energy, dx_filter, dy_filter, dz_filter
		real(prcn) :: kr1, kr2, xi1, eta1, xi2, eta2

		real(prcn), dimension(ndim,ndim) :: tmp_tens, tr1, tr2, bij1, bij2
		real(prcn), allocatable :: resid_stress(:,:,:,:,:)



		real(prcn) :: fluid_frac(nx,my,mz)
		character*50 filename1, filename2
		character*3 tmp_string, string
		integer :: size
		integer :: nx_new, my_new, mz_new

		Real(prcn),  DIMENSION(:,:,:,:), pointer :: ubc_filter, u1uj, u2uj, u3uj

		call screen_separator(30,'^')
		write (*,*) "IN FILTERING"

		if (.not.post_no_flow_mem_alloc) then
			if(.not.allocated(omega)) allocate(omega(nx,my,mz,ndim))
			if(.not.allocated(fr))    allocate(fr(nx,my,mz,ndim))
			if(.not.allocated(ppr))   allocate(ppr(nx,my,mz,ndim))
			if(.not.allocated(diffn)) allocate(diffn(nx,my,mz,ndim))

			ubc_filter => omega
			u1uj       => fr
			u2uj       => ppr
			u3uj       => diffn


			ubc_filter = zero
			u1uj       = zero
			u2uj       = zero
			u3uj       = zero


			call compute_mix_mean_vel
			meanslip_energy = half*mix_mean_slip_mag**2


			filename1 = trim(run_name)//"_dns_field.dat"
			open (unit=1, file=trim(filename1), status="replace")
			write (1,"(2(1A,1I),1A)") "zone I=", nx, " J=", my, " F=point"
			do k=1, mz
				do i=1, nx
					tmp_scal1 = zero
					do dim1=1, ndim
						do dim2=1, ndim
							tmp_tens(dim1,dim2) = (ubc(i,my/2,k,dim1)-ufmean(dim1)) * (ubc(i,my/2,k,dim2)-ufmean(dim2))
						enddo
					enddo
					tmp_tens = tmp_tens / meanslip_energy

					tmp_scal1 = zero
					do dim1=1, ndim
						tmp_scal1 = tmp_scal1 + tmp_tens(dim1,dim1) / 2 
					enddo
					if (fluid_atijk(i,my/2,k)) then
						write (1,"(15d15.7)") i*dx, k*dz, ubc(i,my/2,k,:)/mix_mean_slip_mag, tmp_tens(:,:), tmp_scal1
					else
						write (1,"(15d15.7)") i*dx, k*dz, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
					endif
				enddo
			enddo
			close (1)




			deallocate(nlbc, onlbc)
			allocate(resid_stress(nx,my,mz,ndim,ndim))
			resid_stress = zero

!			filter_length = zero
			fluid_frac    = zero
!			allocate(resid1(nidelta, ndim, ndim), resid2(ndelta, ndim, ndim))


			filename1 = trim(run_name)//"_filtered_values.dat"
			open (unit=1, file=trim(filename1), status="replace")
			close (1)


			do dim1=1, ndim
				do dim2=1, ndim
					do k=1, mz
						do j=1, my
							do i=1, nx
								if (dim1 == 1) then
									u1uj(i,j,k,dim2) = ubcp(i,j,k,dim1) * ubcp(i,j,k,dim2)
								elseif (dim1 == 2) then
									u2uj(i,j,k,dim2) = ubcp(i,j,k,dim1) * ubcp(i,j,k,dim2)
								elseif (dim1 == 3) then
									u3uj(i,j,k,dim2) = ubcp(i,j,k,dim1) * ubcp(i,j,k,dim2)
								endif
							enddo
						enddo
					enddo
				enddo
			enddo
			
		
			do idelta=1, nx, 2
				write (*,*) "DELTA = ", idelta
!				read  (*,*)

				nx_new = nx / idelta
				my_new = my / idelta
				mz_new = mz / idelta

				dx_filter = idelta*dx
				dy_filter = idelta*dy
				dz_filter = idelta*dz

				ubc_filter   = zero
				fluid_frac   = zero
				resid_stress = zero
				do k=1, mz
					do j=1, my
						do i=1, nx

!write (*,*) i,j,k,1

							ib = i-idelta/2
							ie = i+idelta/2

							jb = j-idelta/2
							je = j+idelta/2

							kb = k-idelta/2
							ke = k+idelta/2

							tmp_scal1 = zero
							tmp_vec   = zero
							tmp_tens  = zero

							do kk=kb, ke
								if (kk<1) then
									kkk = kk + mz
								elseif (kk>mz) then
									kkk = kk - mz
								else
									kkk = kk
								endif

								do jj=jb, je
									if (jj<1) then
										jjj = jj + my
									elseif (jj>my) then
										jjj = jj - my
									else
										jjj = jj
									endif

									do ii=ib, ie
										if (ii<1) then
											iii = ii + nx
										elseif (ii>nx) then
											iii = ii - nx
										else
											iii = ii
										endif



!write (*,*) iii,jjj,kkk,2

										if (fluid_atijk(iii,jjj,kkk)) then
											tmp_vec(:) = tmp_vec(:) + ubc(iii,jjj,kkk,:)
											tmp_scal1  = tmp_scal1  + one

											tmp_tens(1,:) = tmp_tens(1,:) + u1uj(iii,jjj,kkk,:)
											tmp_tens(2,:) = tmp_tens(2,:) + u2uj(iii,jjj,kkk,:)
											tmp_tens(3,:) = tmp_tens(3,:) + u3uj(iii,jjj,kkk,:)
										endif
									enddo
								enddo
							enddo

							ubc_filter(i,j,k,:) = tmp_vec(:) / idelta**3
							fluid_frac(i,j,k)   = tmp_scal1  / idelta**3

							tmp_tens(1,:) = tmp_tens(1,:) / idelta**3
							tmp_tens(2,:) = tmp_tens(2,:) / idelta**3
							tmp_tens(3,:) = tmp_tens(3,:) / idelta**3

							do dim1=1, ndim
								do dim2=1, ndim
									resid_stress(i,j,k,dim1,dim2)  = tmp_tens(dim1,dim2) - ubc_filter(i,j,k,dim1) * ubc_filter(i,j,k,dim2)
								enddo
							enddo
						enddo
					enddo
				enddo

				tr1 = zero
				tr2 = zero
				bij1 = zero
				bij2 = zero
				xi1 = zero
				xi2 = zero
				eta1 = zero
				eta2 = zero
				call reynolds_stress_filtered(ubc_filter, resid_stress, tr1, tr2, kr1, kr2, bij1, bij2, xi1, eta1, xi2, eta2)
				resid_stress = resid_stress / meanslip_energy
				tr1 = tr1 / meanslip_energy
				tr2 = tr2 / meanslip_energy
				kr1 = kr1 / meanslip_energy
				kr2 = kr2 / meanslip_energy

				filename1 = trim(run_name)//"_filtered_values.dat"
				open (unit=1, file=trim(filename1), status="old", position="append")
				write (1,"(7f8.4, 42D15.7)") dbydx, lybyd, maxvolfrac, re, rhos/rhof, coeff_rest, idelta*1.0/dbydx,			&
									& tr1(:,:), tr2(:,:), kr1, kr2, &
									& bij1(:,:), bij2(:,:), xi1, eta1, xi2, eta2
				close (1)


				call to_string(idelta,string,size)
				tmp_string = ""
				do i=1, 3-size
					tmp_string="0"//trim(tmp_string)
				enddo
				tmp_string= trim(tmp_string)//trim(string)
				filename2 = trim(run_name)//"_"
				filename2 = trim(filename2)//trim(tmp_string)
				filename2 = trim(filename2)//"_filterd_field.dat"

				open (unit=2, file=trim(filename2), status="replace")
				write (2,"(2(1A,1I),1A)") "zone I=", nx, " J=", my, " F=point"
				do k=1, mz
					do i=1, nx
						tmp_scal1 = zero
						do dim1=1, ndim
							tmp_scal1 = tmp_scal1 + resid_stress(i,my/2,k,dim1,dim1)/2
						enddo
						write (2,"(15d15.7)") i*dx, k*dz, ubc_filter(i,my/2,k,:)/mix_mean_slip_mag, resid_stress(i,my/2,k,:,:), tmp_scal1
					enddo
				enddo
				close (2)
			enddo
			close (1)

			deallocate(resid_stress)
		else
		endif
	end subroutine filter


	subroutine reynolds_stress_filtered(ubc_filter, resid_stress, tr1, tr2, kr1, kr2, bij1, bij2, xi1, eta1, xi2, eta2)
		implicit none
		real(prcn), intent(in)  :: resid_stress(nx,my,mz,ndim,ndim), ubc_filter(nx,my,mz,ndim)
		real(prcn), intent(out) :: tr1(ndim,ndim), tr2(ndim,ndim), kr1, kr2, bij1(ndim,ndim), bij2(ndim,ndim), xi1, eta1, xi2, eta2
		real(prcn) :: ubc_mean(ndim)

		integer :: dim1, dim2
		integer :: i, j, k

		ubc_mean = zero
		do k=1, mz
			do j=1, my
				do i=1, nx
					ubc_mean(:) = ubc_mean(:) + ubc_filter(i,j,k,:)
				enddo
			enddo
		enddo
		ubc_mean(:) = ubc_mean(:) / (nx*my*mz)

		write (*,"(1A,3d15.7)") "THE MEAN SLIP VELOCITY AFTER FILTERATION = ", ubc_mean(:)/mix_mean_slip_mag

		tr1 = zero
		tr2 = zero

		do k=1, mz
			do j=1, my
				do i=1, nx
					do dim1=1, ndim
						do dim2=1, ndim
							tr1(dim1,dim2) = tr1(dim1,dim2) + (ubc_filter(i,j,k,dim1)-ubc_mean(dim1)) * (ubc_filter(i,j,k,dim2)-ubc_mean(dim2))
						enddo
					enddo
				enddo
			enddo
		enddo

		do k=1, mz
			do j=1, my
				do i=1, nx
!					if (fluid_atijk(i,j,k)) then
!						tr1(:,:) = tr1(:,:) + resid_stress(i,j,k,:,:)
!					endif
					tr2(:,:) = tr2(:,:) + resid_stress(i,j,k,:,:)
				enddo
			enddo
		enddo

!		tr1(:,:) = tr1(:,:) / count_fluid
		tr1(:,:) = tr1(:,:) / (nx*my*mz)
		tr2(:,:) = tr2(:,:) / (nx*my*mz)

		kr1 = zero
		kr2 = zero
		do i=1, ndim
			kr1 = kr1 + tr1(i,i)/2
			kr2 = kr2 + tr2(i,i)/2
		enddo

		bij1 = zero
		bij2 = zero
		call calc_anisotropy(tr1, bij1, xi1, eta1)
		call calc_anisotropy(tr2, bij2, xi2, eta2)
	end subroutine reynolds_stress_filtered

	SUBROUTINE calc_anisotropy(uij, bij, xi, eta)
		IMPLICIT NONE

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










  SUBROUTINE flow_snapshot4
    Use nlarrays , only : ur1,uf1
    USE dem_mod, only : is_mobile, des_pos_new, des_radius
    IMPLICIT NONE
    Integer  :: sunit,i,j,k,l,m,isp, mark, idim
    INTEGER, SAVE :: zone_count = 0
    LOGICAL, SAVE :: first_pass=.TRUE.
    REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
         & mean_vel(ndim), mean_force(ndim), position(ndim)
    CHARaCTER*50 :: FILENAME1, filename2, filename3
    integer, save :: sphrunit1, sphrunit2, sphrunit3
    CHARACTER(LEN=80) :: formfile

	real(prcn), allocatable :: out_arr(:,:,:,:), trans_buf(:)
	integer :: node_num, iproc
	integer :: jj, j1, j2, j3
	real(prcn) :: tmp, tmp1, tmp2, tmp3
	real(prcn), pointer, dimension(:,:,:,:) :: ur
	

	call screen_separator(30,'^')
	write (*,*) "IN SNAPSHOT...."

    formfile='formatted' 
    sunit  = getnewunit(minunitno,maxunitno)


	ur=>ubcp


!	j1=1
!	j2=my/2
!	j3=my
!

!	allocate(out_arr(mx1,3,mz,ndim+1))
!	do idim=1, ndim
!		do k=1, mz
!			do jj=1, 3
!				if (jj==1) then
!					j=j1
!				elseif (jj==2) then
!					j = j2
!				elseif (jj==3) then
!					j = j3
!				endif
!
!				do i=1, mx1
!					out_arr(i,jj,k,idim) = velr(i,j,k,idim)-ufmean(idim)
!!					if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
!				enddo
!			enddo
!		enddo
!	enddo

	if (I_AM_NODE_ZERO) then
		write (*,*) "GENERATING THE SNAPSHOT OF THE FIELD"
!		out_arr(:,:,:,:) = out_arr(:,:,:,:) / umeanslip
!
!		do k=1, mz
!			do j=1, 3
!				do i=1, mx1
!					tmp = 0d0
!					do idim=1, ndim
!						tmp = tmp + out_arr(i,j,k,idim) * out_arr(i,j,k,idim)
!					enddo
!					tmp = tmp/2
!					out_arr(i,j,k,4) = tmp
!				enddo
!			enddo
!		enddo

		sphrunit1 = 31 !getnewunit(minunitno,maxunitno)
!		sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
!		sphrunit3 = 33 !getnewunit(minunitno,maxunitno)

		FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion.dat'
!		FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'
!		filename2 = TRIM(RUN_NAME)//'_sphr_motion_act1.dat'
!		filename3 = TRIM(RUN_NAME)//'_sphr_motion_act2.dat'


		IF(first_pass)THEN
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='replace')
			first_pass = .FALSE.

!			sphrunit  = 31 !getnewunit(minunitno,maxunitno)
!			sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
!			sphrunit3 = 33 !getnewunit(minunitno,maxunitno)
!
!			FILENAME = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'
!			filename2 = TRIM(RUN_NAME)//'_sphr_motion_act1.dat'
!			filename3 = TRIM(RUN_NAME)//'_sphr_motion_act2.dat'

!			CALL  RUN_TIME_FILE_OPENER(sphrunit,FILENAME, formfile)
!			CALL  RUN_TIME_FILE_OPENER(sphrunit2,FILENAME2, formfile)
!			CALL  RUN_TIME_FILE_OPENER(sphrunit3,FILENAME3, formfile)

			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='replace')
!			OPEN(unit = sphrunit2,file=TRIM(FILENAME2),form='formatted', status='unknown')
!			OPEN(unit = sphrunit3,file=TRIM(FILENAME3),form='formatted', status='unknown')
		ELSE
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', POSITION='append')


			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='old', position='append')
!			OPEN(unit = sphrunit2,file=TRIM(FILENAME2),form='formatted', status='old', position='append')
!			OPEN(unit = sphrunit3,file=TRIM(FILENAME3),form='formatted', status='old', position='append')
		endif

		write(sunit,*)'ZONE T = "', t, '",'
		write(sunit,*)'DATAPACKING=POINT, I =', mx1,  ', J=', my, ', K=', mz


		tmp1 = half*dot_product(ufmean(:)-usmean(:), ufmean(:)-usmean(:))


		do k=1, mz
			do j=1, my
				do i=1, mx1
					tmp2 = dot_product(ur(i,j,k,:)-ufmean(:), ur(i,j,k,:)-ufmean(:))/2 / tmp1
					write(sunit,*) i, j, k, tmp2 !out_arr(i,jj,k,1), out_arr(i,jj,k,2), out_arr(i,jj,k,3), out_arr(i,jj,k,4)
				enddo
			enddo
		enddo
		close(sunit,status='keep')

		WRITE(sphrunit1,*)'ZONE T= "', t, ' " '
!		WRITE(sphrunit2,*)'ZONE T= "', t, ' " '
!		WRITE(sphrunit3,*)'ZONE T= "', t, ' " '


!		close(sphrunit1)
!		close(sphrunit2)
!		close(sphrunit3)
!		stop

		do idim = 1, ndim
			mean_force(idim) = SUM(force(1:nbody,idim))/real(nbody,prcn)
			mean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
		end do

		DO m=1,nbody
			fluct_vel(1:ndim) = velbdy(m,1:ndim)-mean_vel(1:ndim)
			fluct_force(1:ndim) = force(m,1:ndim)-mean_force(1:ndim)
			mark = -1
			do idim = 1, ndim
				position(idim) = XC(m,idim) !+frame_pos(idim)
				if(position(idim).lt.one) position(idim) = position(idim) + real(my,prcn)
				if(position(idim).ge.real(my+1,prcn))position(idim) = position(idim) - real(my,prcn)
			end do

!			tmp1 = j1*dy
!			tmp2 = j2*dy
!			tmp3 = j3*dy

!			if (abs(xc(m,2)-j1)<=radbdy(m).or.abs(xc(m,2)-j3)<=radbdy(m)) then
!				WRITE(sphrunit2,'(10(2x,f12.8))')  position(1), position(2), position(3), radbdy(m),radbdy(m),real(mark)
!			elseif (abs(xc(m,2)-j2)<=radbdy(m)) then
!				WRITE(sphrunit3,'(10(2x,f12.8))')  position(1), position(2), position(3), radbdy(m),radbdy(m),real(mark)
!			else
				WRITE(sphrunit1,'(10(2x,f12.8))')  position(1), position(2), position(3), radbdy(m) !,radbdy(m),real(mark)
!			endif
		enddo
		close (sphrunit1)
!		close (sphrunit2)
!		close (sphrunit3)
!		deallocate(out_arr)
	endif
	call screen_separator(30,'-')
  END SUBROUTINE flow_snapshot4



	subroutine c_k
		implicit none
		real(prcn) :: rhom, tp, tf1, tf2, tf3, stokes1, stokes2, stokes3, phi, es, ef, em, ck1, ck2, ck3, ck4
		integer :: nvar, nvar1, nvar2
		character*50 filename

		if (.not.post_no_flow_mem_alloc) then
			write (*,*) "C_K..."
			rhom = rhof*(1-maxvolfrac) + rhos*maxvolfrac
			tp  = rhos/rhof * char_length**2 / 18/vis
			tf1 = char_length / mix_mean_slip_mag
			tf2 = sqrt(vis/sijsij_pure)
			tf3 = tke_pure / sijsij_pure

			stokes1 = tp / tf1
			stokes2 = tp / tf2
			stokes3 = tp / tf3

			ef = rhof*(1-maxvolfrac)*tke_pure
			es = rhos* (maxvolfrac) *tke_s_pure(1)
			em = ef + es
			phi = rhos*maxvolfrac / (rhof*(1-maxvolfrac))

			ck1 = phi * (tke_s_pure(1)/tke_pure) / (1+phi*(tke_s_pure(1)/tke_pure))
			ck2 = phi / (1+phi+stokes1)
			ck3 = phi / (1+phi+stokes2)
			ck4 = phi / (1+phi+stokes3)

			if (I_AM_NODE_ZERO) then
				filename = trim(run_name)//"_ck.dat"
				open (unit=1, file=trim(filename), status="replace")

				if (.not.imove==1) then
					if (nphases==1) then
						write (1,"(4f8.4,12D15.7)") dbydx, lybyd, maxvolfrac, re,								&
											&	phi, tp, tf1, tf2, stokes1, stokes2, ef, es, em, ck1, ck2, ck3
					endif
				else
					if (nphases==1) then
						write (1,"(7f8.4,18D15.7)") dbydx, lybyd, maxvolfrac, re, rhos/rhof, rhom/rhof, coeff_rest,	&
											&	phi, tke_pure, tke_s_pure(1), tke_pure/tke_s_pure(1), tp, tf1, tf2, tf3, &
											&	stokes1, stokes2, stokes3, ef, es, em, ck1, ck2, ck3, ck4
					endif
				endif
				close (1)
			endif
		else
			if(I_AM_NODE_ZERO)then
				if (.not.imove==1) then
					if (nphases==1) then
						nvar  = 19
						nvar1 = 4
						nvar2 = 15
					endif
				else
					if (nphases==1) then
						nvar  = 25
						nvar1 = 7
						nvar2 = 18
					endif
				endif
				filename = "_ck.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line)
			endif 
		endif	
	end subroutine c_k

	subroutine pi_groups
		implicit none
		integer :: iphs, m, idim, partstart, partend
		real(prcn) :: grant(nphases), stokes, pres_part, visc_part, int_tke, tmp1, pres_part_fluc, visc_part_fluc, int_tke_fluc, &
				&	gamma_coll, gamma_vis, gamma_rat1, gamma_rat2, kay, r_diss, force(ndim), force_mag, tfluid, tcoll, g_r, g

		real(prcn) :: p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, tp, tf, tt, tt1, tt2, int_tke_tot
		integer :: nvar, nvar1, nvar2
		character*50 filename


		if (.not.post_no_flow_mem_alloc) then
#if 0
			write (*,*) "IN PI_GROUPS"
			if (.not.allocated(acc_mean)) then
				write (*,*) "FIRST RUN AiVi ..."
				stop
			endif

			grant = zero
			do iphs = 1, nphases
				partstart = phase_array(iphs)%pstart
				partend = phase_array(iphs)%pend

				do m = partstart, partend
					do idim = 1, ndim
						grant(iphs) = grant(iphs) + (phase_array(iphs)%mean_spec_vel(idim)-velbdy(m,idim))**2.d0
					end do
				end do
				grant(iphs) = grant(iphs)/(three*phase_array(iphs)%npart)
			end do

			if (nphases==1) then
				ReT = char_length*sqrt(grant(1))/vis
				write (*,"(1A,2D15.7)") "T and ReT = ", grant(1), ReT
			endif

			stokes = rhos/rhof * ret/9
			write (*,"(1A,1D15.7)") "St(ReT) = ", stokes


			tt1 = sqrt(3*grant(1)) / acc_mean(1)
			tt2 = sqrt(3*grant(1)) / acc_var(1)

			if (abs(phiavg-0.1)<0.01) then
				g_r = 1.296
			elseif (abs(phiavg-0.2)<0.01) then
				g_r = 1.763
			elseif (abs(phiavg-0.3)<0.01) then
				g_r = 2.329
			elseif (abs(phiavg-0.4)<0.01) then
				g_r = 3.343
			endif

			write (*,"(1A,1D15.7)") "g(r) = ", g_r

			g = 8 * char_length * nbody / (doml(1)*doml(2)*doml(3)) * sqrt(6. / (rhos*char_length**3)) * g_r
			tcoll = 2. / (g*sqrt(grant(1)))

			write (*,"(1A,2d15.7)") "<A>, s_A, T_vis = ", acc_mean(1), acc_var(1)
			write (*,"(1A,3d15.7)") "T_<A>, T_sA, T_coll = ", tt1, tt2, tcoll

			filename = trim(run_name)//"_scales.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")

			write (1,"(6f8.4,10D15.7)") dbydx, lybyd, maxvolfrac, re, rhos, coeff_rest,	&
				&	ReT, grant(1)/mix_mean_slip_mag**2, rhos/rhof*3*grant(1)/mix_mean_slip_mag**2, acc_mean(1), acc_var(1), tt1*mix_mean_slip_mag/char_length, &
				&	tt2*mix_mean_slip_mag/char_length, tcoll*mix_mean_slip_mag/char_length, tt1/tcoll, tt2/tcoll

			close (1)
		else
			nvar  = 16
			nvar1 = 6
			nvar2 = 10

			filename = "_scales.dat"
			call mis_average(nvar, nvar1, nvar2, filename, line)
		endif
#endif



#if 0
			p1 = maxvolfrac
			p2 = coeff_rest
			p3 = rhos/rhof
			p4 = sqrt(tke_pure)/mix_mean_slip_mag
			p5 = sqrt(grant(1))/mix_mean_slip_mag
			p6 = mix_mean_slip_mag * char_length / vis

			tp = rhos/rhof * char_length**2 / 18/vis
			tf = char_length / mix_mean_slip_mag
			p7 = tp / tf

			tt1 = sqrt(3*grant(1)) / acc_mean(1)
			tt2 = sqrt(3*grant(1)) / acc_var(1)

			if (abs(phiavg-0.1)<0.01) then
				g_r = 1.296
			elseif (abs(phiavg-0.2)<0.01) then
				g_r = 1.763
			elseif (abs(phiavg-0.3)<0.01) then
				g_r = 2.329
			elseif (abs(phiavg-0.4)<0.01) then
				g_r = 3.343
			endif

			write (*,"(1A,1D15.7)") "g(r) = ", g_r

			g = 8 * char_length * nbody / (doml(1)*doml(2)*doml(3)) * sqrt(6. / (rhos*char_length**3)) * g_r
			tcoll = 2. / (g*sqrt(grant(1)))
			p8 = tt1 / tcoll

			write (*,"(1A,3d15.7)") "SQRT(3T), <A>, T_vis = ", sqrt(3*grant(1)), acc_mean(1), tt
			write (*,"(1A,3d15.7)") "T_vis, T_coll, RATIO = ", tt, tcoll, p8

			gamma_coll = 24./(char_length*sqrt(pi)) * (1.-coeff_rest) * rhos * maxvolfrac**2 * g_r * grant(1)**1.5

			p9 = gamma_coll /rhos /(grant(1)/tcoll)

			gamma_vis = dissipation
			p10 = gamma_vis /(grant(1)/tcoll)

			p11 = sijsij_pure / (vis*(mix_mean_slip_mag/char_length)**2)

			p12 = norm_drag_spec(1)

			gamma_rat1 = gamma_coll/rhos/gamma_vis
			gamma_rat2 = gamma_coll/rhos/sijsij_pure
!			gamma_rat2 = p9/p10 !* p5**2 * p6/p8
#endif

#if 0
			pres_part = zero
			visc_part = zero
			int_tke   = zero
			force = zero
			do m=1, nbody
				do idim=1, ndim
					tmp1 = ufmean(idim)-usmean(idim)

					pres_part = pres_part + tmp1 * pres(m,idim)
					visc_part = visc_part + tmp1 * visc(m,idim)
					int_tke   = int_tke   + tmp1 * (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)

					force(idim) = force(idim) + (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)
				enddo
			enddo
			force_mag = sqrt(dot_product(force,force))

			pres_part_fluc = zero
			visc_part_fluc = zero
			int_tke_fluc   = zero
			do m=1, nbody
				do idim=1, ndim
					tmp1 = velbdy(m,idim) - usmean(idim)

					pres_part_fluc = pres_part_fluc - pres(m,idim)                * tmp1
					visc_part_fluc = visc_part_fluc - visc(m,idim)                * tmp1
					int_tke_fluc   = int_tke_fluc   - (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0) * tmp1
				enddo
			enddo

			pres_part = pres_part / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			visc_part = visc_part / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			int_tke   = int_tke   / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac))

			pres_part_fluc = pres_part_fluc / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			visc_part_fluc = visc_part_fluc / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			int_tke_fluc   = int_tke_fluc   / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac))

			int_tke_tot = (int_tke + int_tke_fluc)
#endif


			write (*,"(4d15.7)") gamma_coll / int_tke_tot, gamma_vis*rhos / int_tke_tot, (gamma_coll+gamma_vis*rhos) / int_tke_tot, gamma_rat1

!			write (*,"(1A,1d15.7)") "GAMMA_COLL/GAMMA_VIS = ", gamma_rat1
!			write (*,"(1A,1d15.7)") "GAMMA_COLL/DISIP_VIS = ", gamma_rat2
!			write (*,"(1A,1d15.7)") "GAMMA_(COLL+VIS)/DISIP_VIS = ", (gamma_coll+rhos*gamma_vis) / (int_tke + int_tke_fluc)

			write (*,"(1A,1d15.7)") "TIME_RAT = ", p8

			filename = trim(run_name)//"_pi_groups.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")

			write (1,"(6f8.4,15D15.7)") dbydx, lybyd, maxvolfrac, re, rhos, coeff_rest,	&
				&	p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ReT, gamma_rat1, gamma_rat2

			close (1)


			filename = trim(run_name)//"_energy_balance.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")
			write (1,"(6f8.4,6D15.7)") dbydx, lybyd, maxvolfrac, re, rhos, coeff_rest,	&
				&	int_tke_pure, sijsij_pure, gamma_coll, &
				&	int_tke_pure/ (vis*(mix_mean_slip_mag/char_length)**2), sijsij_pure / (vis*(mix_mean_slip_mag/char_length)**2), &
				&	gamma_coll/ (vis*(mix_mean_slip_mag/char_length)**2)
			close (1)
		else
!			nvar  = 21
!			nvar1 = 6
!			nvar2 = 15

!			filename = "_pi_groups.dat"
!			call mis_average(nvar, nvar1, nvar2, filename, line)

			nvar  = 12
			nvar1 = 6
			nvar2 = 6

			filename = "_energy_balance.dat"
			call mis_average(nvar, nvar1, nvar2, filename, line)
		endif

#if 0
			pres_part = zero
			visc_part = zero
			int_tke   = zero
			force = zero
			do m=1, nbody
				do idim=1, ndim
					tmp1 = ufmean(idim)-usmean(idim)

					pres_part = pres_part + tmp1 * pres(m,idim)
					visc_part = visc_part + tmp1 * visc(m,idim)
					int_tke   = int_tke   + tmp1 * (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)

					force(idim) = force(idim) + (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)
				enddo
			enddo
			force_mag = sqrt(dot_product(force,force))

			pres_part_fluc = zero
			visc_part_fluc = zero
			int_tke_fluc   = zero
			do m=1, nbody
				do idim=1, ndim
					tmp1 = velbdy(m,idim) - usmean(idim)

					pres_part_fluc = pres_part_fluc - pres(m,idim)                * tmp1
					visc_part_fluc = visc_part_fluc - visc(m,idim)                * tmp1
					int_tke_fluc   = int_tke_fluc   - (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0) * tmp1
				enddo
			enddo

			pres_part = pres_part / nbody !(doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			visc_part = visc_part / nbody !(doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			int_tke   = int_tke   / nbody !(doml(1)*doml(2)*doml(3)*(1-maxvolfrac))

			pres_part_fluc = pres_part_fluc / nbody !(doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			visc_part_fluc = visc_part_fluc / nbody !(doml(1)*doml(2)*doml(3)*(1-maxvolfrac))
			int_tke_fluc   = int_tke_fluc   / nbody !(doml(1)*doml(2)*doml(3)*(1-maxvolfrac))


			r_diss = (int_tke + int_tke_fluc) / (9*pi*vis*char_length*grant(1))

			kay = (1-maxvolfrac/2)/(1-maxvolfrac)**3
			gamma_rat = rhos/rhof * 4 * (1-coeff_rest) * maxvolfrac * kay * ReT / (9*sqrt(pi)*r_diss)

			write (*,"(1A,2d15.7)") "TKE_TRANSFER, MEAN & FLUC = ", int_tke, int_tke_fluc
			write (*,"(1A,1d15.7)") "R_DISS = ", r_diss

			write (*,"(1A,2d15.7)") "KAY = ", kay
			write (*,"(1A,2d15.7)") "GAMMA_RAT = ", gamma_rat



			tcoll  = char_length/(24*maxvolfrac*kay)*sqrt(pi/grant(1))
			tfluid = pi*((two*radbdy(1)*dx)**3.d0)/6.d0*rhos *sqrt(3*grant(1))/force_mag

			write (*,"(1A,1d15.7)") "TIME_RAT = ", tfluid/tcoll


			filename = trim(run_name)//"_dissip_rat.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")

			write (1,"(6f8.4,3D15.7)") dbydx, lybyd, maxvolfrac, re, rhos, coeff_rest,	ReT, tfluid/tcoll, gamma_rat
			close (1)
		else
			nvar  = 9
			nvar1 = 6
			nvar2 = 3

			filename = "_dissip_rat.dat"
			call mis_average(nvar, nvar1, nvar2, filename, line)
		endif
#endif

	end subroutine pi_groups


	subroutine neighb_dist
		implicit none
		integer :: i, j, idim, ibin, nearbin, nnear
		integer :: nbin=200
		real(prcn), allocatable :: h_near(:), h_t(:)
		real(prcn) :: rmax, rad(ndim), radmod, rad_reg, rbin, dr, total_n, total_t
		logical :: nearest
		character*50 filename


		rmax = lybyd*half
		dr = rmax / nbin

		allocate(h_near(nbin), h_t(nbin))
		h_near = zero
		h_t    = zero

		call screen_separator(30,'^')
		write (*,*) "IN NEIGHB_DIST"

!		write (*,*) lybyd, rmax, dr, nbody
!		read (*,*) 

		write (*,*) "RAD_REG = ", radbdy(1)*2*rad_factor*dx

		do i=1, nbody-1
			rad_reg = radbdy(i)*2*rad_factor*dx

			nearbin = 1000
			nnear = 0

			do j=i+1, nbody
				rad(:) = abs(xc(i,:)-xc(j,:)) * dx

				do idim=1, 3
					if (rad(idim) > rmax) rad(idim) = lybyd - rad(idim)
				enddo
				radmod = sqrt(rad(1)**2 + rad(2)**2 + rad(3)**2)

				if (radmod <= rmax) then
					ibin = int(radmod / dr) + 1
					rbin = (ibin-0.5) * dr

					if (ibin == nearbin) then
						nnear = nnear + 1
					elseif (ibin < nearbin) then
						nearbin = ibin
						nnear = 1
					endif

					if (rbin <= rad_reg) then
						h_t(ibin) = h_t(ibin) + 1
					endif
				endif
			enddo

			if (nearbin<nbin) h_near(nearbin) = h_near(nearbin) + nnear
		enddo

!		h(:) = h(:) / nbody
		total_n = sum(h_near(:))
		total_t = sum(h_t(:))
		write (*,*) "TOTAL SUM H_NEAR = ", total_n
		write (*,*) "NORMALIZING WITH = ", total_n*dr
		write (*,*)
		write (*,*) "TOTAL SUM H_TOT  = ", total_t
		write (*,*) "NORMALIZING WITH = ", total_t*dr
		write (*,*)
		write (*,*) "WRITING OUTPUT"
		call screen_separator(30,'-')

		filename = trim(run_name)//"_nearest_h.dat"
		open (unit=1,file=trim(filename),status="replace",action="write")		

		do i=1, nbin
			write (1,"(3D15.7, 2F6.2)") (i-0.5) * dr, h_near(i)/total_n/dr, h_t(i)/total_t/dr, lybyd, maxvolfrac
		enddo
		close (1)

		deallocate(h_near, h_t)

	end subroutine neighb_dist



	subroutine dissip_continuous
		implicit none

		character*50 filename
		integer :: m, n, l, ith, iphs


		real(prcn) :: x,y,z,rad,tmp1, th, dth
		integer :: nth, nbnd

		real(prcn), allocatable :: dissip_visc(:), dissip_pres1(:), dissip_pres2(:)
		integer, allocatable :: count(:)



		integer, dimension(ndim) :: vcellb !, vcello, vcello2
		integer, dimension(ndim) :: pcellb !, pcello, pcello2
		real(prcn), dimension(ndim) :: xl 	!, xlo, xlo2
		real(prcn), dimension(ndim) :: xpb 	!, xlop, xlo2p
		real(prcn), dimension(ndim) :: vel, ul 	!, ulo, ulo2, ul_p, ulo_p, ulo2_p
		real(prcn), dimension(ndim) ::nll, onll ,dfll, pgl, pgl1, pgl2	!, pglo, pglo2
		real(prcn) :: pl					!, plo, plo2


		integer :: ib, ie, jb, je, kb, ke, onew




		write (*,*) "COMPUTING THE DISSIPATION"

		tmp1 = 0.5*phase_array(nphases)%dia/doml(2)*my
		nth  = nint(twopi*tmp1*f2/4)*2

		write (*,*) "NUMBER OF BINS: ", nth

!		call calc_visc
!		call calc_pgrad


		do m=1, nbody
			iphs = 1
			nbnd = phase_array(iphs)%nbnd
			NULLIFY(bndarray)
			bndarray => phase_array(iphs)%bndpts

			tmp1 = 0.5*phase_array(nphases)%dia/doml(2)*my
			nth  = nint(twopi*tmp1*f2/4)*2

			allocate(dissip_visc(0:nth), dissip_pres1(0:nth), dissip_pres2(0:nth), count(0:nth))
			dissip_visc = zero
			dissip_pres1 = zero
			dissip_pres2 = zero
			count  = 0

			DO l=1,nbnd
				x = bndarray(1,l)
				y = bndarray(2,l)
				z = bndarray(3,l)

				rad = sqrt(x**2+y**2+z**2)

				th  = acos(abs(x)/rad)
				if (x>0) th = pi-th

				dth = twopi/2/nth
				ith = nint(th/dth)
				th  = dth*ith

				if (th<0.or.th>pi) write (*,*) "ERROR IN THETA"
				if (ith>nth)       write (*,*) "ERROR IN NTHETA"

!				nphi=nint(twopi*tmp1*sin(th)*f2/4)*4
!				if (abs(sin(th))>small_number) then
!					phi = acos(abs(y/rad/sin(th)))
!					if (abs( y / rad / sin(th)) > one) phi=zero
!
!					if ( y >= zero .and. z >= zero) then
!
!					elseif (y < zero .and. z >= zero) then
!						phi = pi - phi
!					elseif ( y < zero .and. z < zero ) then
!						phi = pi + phi
!					elseif ( y >= zero .and. z < zero) then
!						phi = twopi - phi
!					endif
!
!					dphi = twopi/nphi
!					iphi = nint(phi/dphi)
!!					if (iphi>=nphi) iphi = 0
!					phi  = dphi*iphi
!				else
!					phi  = 0
!					nphi = 1
!					iphi = 0
!				endif

				rad = zero
				DO n=1,ndim
					xl(n) = xc(m,n)+ bndarray(n,l)*radbdy(m)
					rad   = rad+(bndarray(n,l)*radbdy(m))**2.0
!					is(n) = INT(xl(n))
					ul(n) = zero
					pgl(n)= zero
				enddo
				rad = DSQRT(rad)
				xpb(1) = xl(1)-0.5
				xpb(2:3)=xl(2:3)
				do n = 1, ndim
					if(xpb(n).lt.zero) then 
						pcellb(n) = int(xpb(n)-1)
					else 
						pcellb(n) = int(xpb(n))
					end if
					if(xl(n).lt.zero) then 
						vcellb(n) = int(xl(n)-1)
					else 
						vcellb(n) = int(xl(n))
					end if
				end do
#if 0
#if PARALLEL
				xltemp = xl(1)
				xptemp = xpb(1)
				pcelltemp = pcellb(1)
				vcelltemp = vcellb(1)
				if(l.eq.FOCUS_POINT)then
					PRINT*,' xl = ', myid, xltemp, xptemp, pcelltemp
				end if
				if(.not.CELL_IN_VEL_GRID(vcelltemp))then
					WEST_PERIODIC_IMAGE(vcellb(1),vcelltemp,xl(1),xltemp)
					WEST_PERIODIC_IMAGE(pcellb(1),pcelltemp,xpb(1),xptemp)
					EAST_PERIODIC_IMAGE( vcellb(1),vcelltemp,xl(1),xltemp)
					EAST_PERIODIC_IMAGE_PRES(pcellb(1),pcelltemp,xpb(1),xptemp)
					if(l.eq.FOCUS_POINT)then
						PRINT*,' xl IMAGES = ', myid,xltemp, xptemp, pcelltemp
					end if

					if(.not.CELL_IN_VEL_GRID(vcelltemp)) goto 777
				end if
				vcellb(1) = vcelltemp
				pcellb(1) = pcelltemp
				xl(1) = xltemp
				xpb(1) = xptemp
#endif
#endif
				pl=zero
				call interpolate_pdata(pcellb,xpb,pgl,pl,l)
				call interpolate_udata(vcellb,xl,ib,ie,jb,je,kb,ke,ul,nll,onll,dfll, 1, m, l, onew)

				vel = usmean-ufmean

				pgl1 = pgl
				pgl2 = pgl - mpg

				dissip_visc(ith) = dissip_visc(ith) + dot_product(vel, dfll)
				dissip_pres1(ith) = dissip_pres1(ith) + dot_product(vel, pgl1)
				dissip_pres2(ith) = dissip_pres2(ith) + dot_product(vel, pgl2)
				count(ith)  = count(ith) + 1
	 		enddo
		enddo

		do ith=1, nth
			if (count(ith) > 0) then
				dissip_visc(ith) = dissip_visc(ith) / count(ith)
				dissip_pres1(ith) = dissip_pres1(ith) / count(ith)
				dissip_pres2(ith) = dissip_pres2(ith) / count(ith)
			endif
		enddo

		dissip_visc  =  dissip_visc / (vis*(umeanslip/dia_phys)**2)
		dissip_pres1 = -dissip_pres1 / (vis*(umeanslip/dia_phys)**2)
		dissip_pres2 = -dissip_pres2 / (vis*(umeanslip/dia_phys)**2)

		write (*,*) vis, umeanslip, dia_phys, (vis*(umeanslip/dia_phys)**2)
		write (*,*) mpg

		write (*,*) "WRITING THE OUTPUT..."

		filename = trim(run_name)//"_dissip_local_cont1.dat"
		open (unit=1,file=trim(filename),status="replace",action="write")		
		do ith=1, nth
			write (1,"(1F6.2, 5D15.7)") dth*ith*180/pi, dissip_visc(ith), dissip_pres1(ith), dissip_pres2(ith), dissip_visc(ith)+dissip_pres1(ith), dissip_visc(ith)+dissip_pres2(ith)
		enddo
		close (1)
	end subroutine dissip_continuous

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

	subroutine compute_AiVi
		implicit none
		real(prcn), allocatable :: acc_fluc(:,:), vel_fluc(:,:), acc_var(:,:), vel_var(:,:), acc_fluc_meanf(:,:), acc_var_meanf(:,:)

		allocate(acc_fluc(nbody,ndim), vel_fluc(nbody,ndim), acc_var(nphases,ndim), vel_var(nphases,ndim))
		allocate(acc_fluc_meanf(nbody,ndim), acc_var_meanf(nphases,ndim))

		call AiVi(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)

		deallocate(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)
	end subroutine compute_AiVi


	subroutine AiVi(acc_fluc, vel_fluc, acc_var, vel_var, acc_fluc_meanf, acc_var_meanf)
		implicit none

		real(prcn), intent(out) :: acc_fluc(nbody,ndim), vel_fluc(nbody,ndim), acc_var(nphases,ndim), vel_var(nphases,ndim), acc_fluc_meanf(nbody,ndim), acc_var_meanf(nphases,ndim)

		integer :: pstart, pend, m, iphs, idim
		character*50 filename
		real(prcn), allocatable :: acc_avg(:,:), acc_avg_meanf(:,:), phase_mass(:)
		real(prcn) :: tmp1, re_tmp, f_tmp
		integer :: nvar1, nvar2, nvar
		logical :: filexist

		if (I_AM_NODE_ZERO) then

		if (.not.post_no_flow_mem_alloc) then
			write (*,*) "IN AiVi..."
			if (.not.allocated(acc_mean)) allocate(acc_mean(nphases), vel_mean(nphases))
			if (.not.allocated(acc_avg_meanf)) allocate(acc_avg_meanf(nphases,ndim))
			if (.not.allocated(acc_avg)) allocate(acc_avg(nphases,ndim), phase_mass(nphases))

			! PARTICLE ACCELERATION AND ITS SD OBTAINED FROM THE MEAN DRAG FORCE MODEL
			acc_avg_meanf = 0
			pstart = 1
			do iphs = 1, nphases
				phase_mass(iphs) = rhos*pi*(phase_array(iphs)%dia)**3.d0/6.d0
				pend = pstart + phase_array(iphs)%npart- 1
				do m = pstart, pend
					re_tmp = sqrt(dot_product(ufmean(:)-velbdy(m,:), ufmean(:)-velbdy(m,:)))
					re_tmp = re_tmp * (1-maxvolfrac) * dia_phys / vis

					call compute_ibm_drag(maxvolfrac, re_tmp, f_tmp)
					acc_fluc_meanf(m,:) = f_tmp * 3.d0 * pi * vis * dia_phys * (1-maxvolfrac) * (ufmean(:)-velbdy(m,:)) / phase_mass(iphs)
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


!			MEAN ACCELERATION
			acc_avg = zero
			pstart = 1
			do iphs = 1, nphases
				phase_mass(iphs) = rhos*pi*(phase_array(iphs)%dia)**3.d0/6.d0
				pend = pstart + phase_array(iphs)%npart- 1
				do m = pstart, pend
					acc_avg(iphs,:) = acc_avg(iphs,:)+force(m,:)/phase_mass(iphs)
				enddo
				acc_avg(iphs,:) = acc_avg(iphs,:) / (pend-pstart+1)
				pstart = pend + 1
			end do

			do iphs = 1, nphases
				acc_avg(iphs,:) = acc_avg(iphs,:) ! + mpg(:)*pi*((two*radbdy(1)*dx)**3.d0)/6.d0/phase_mass(iphs)
			enddo

			do iphs = 1, nphases
				acc_mean(iphs) = sqrt(dot_product(acc_avg(iphs,:), acc_avg(iphs,:)))
				vel_mean(iphs) = sqrt(dot_product(phase_array(iphs)%mean_spec_vel(:), phase_array(iphs)%mean_spec_vel(:)))
			enddo

	!		FLUCTUATING ACCELERATION
			acc_fluc = zero
			pstart = 1
			do iphs = 1, nphases
				pend = pstart + phase_array(iphs)%npart- 1
				do m = pstart, pend
					do idim=1, ndim
						acc_fluc(m,idim) = force(m,idim)/phase_mass(iphs) - acc_avg(iphs,idim)
					end do
				enddo
				pstart = pend + 1
			end do

	!		FLUCTUATING VELOCITY
			vel_fluc = zero
			pstart = 1
			do iphs = 1, nphases
				pend = pstart + phase_array(iphs)%npart- 1
				do m = pstart, pend
					do idim=1, ndim
						vel_fluc(m,idim) = velbdy(m,idim)-phase_array(iphs)%mean_spec_vel(idim)
					end do
				enddo
				pstart = pend + 1
			end do

	!		STANDARD DEVIATIONS
			acc_var = zero
			vel_var = zero

			do idim=1, ndim
				pstart = 1
				do iphs = 1, nphases
					pend = pstart + phase_array(iphs)%npart- 1
					do m = pstart, pend
						acc_var(iphs,idim) = acc_var(iphs,idim) + (force(m,idim)/phase_mass(iphs) - acc_avg(iphs,idim))**2.d0
						vel_var(iphs,idim) = vel_var(iphs,idim) + (velbdy(m,idim) - phase_array(iphs)%mean_spec_vel(idim))**2.d0
					enddo

					acc_var(iphs,idim) = acc_var(iphs,idim) / (pend-pstart+1) !/ ndim
					vel_var(iphs,idim) = vel_var(iphs,idim) / (pend-pstart+1) !/ ndim

					acc_var(iphs,idim) = sqrt(acc_var(iphs,idim))
					vel_var(iphs,idim) = sqrt(vel_var(iphs,idim))

					pstart = pend + 1
				end do
			enddo

!write (*,"(3(3d15.7))") acc_var(:,:)
!write (*,"(3(3d15.7))") vel_var(:,:)


	!		SOURCE / DISSIPATION
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
			source      = two/three * source      / real(nbody, prcn)
			dissipation = two/three * dissipation / real(nbody, prcn)

			if (I_AM_NODE_ZERO) then
				filename=trim(run_name)//"_AiVi"
				if (from_post) then
					filename=trim(filename)//"_post.dat"
					open (unit=1, file=trim(filename), status="replace", action="write")
				else
					filename=trim(filename)//".dat"
					if (irestart==0.and.first_pass) then
						open (unit=1, file=trim(filename), status="replace", action="write")
					else
						inquire (file=trim(filename), exist=filexist)
						if (.not.filexist) then
							open (unit=1, file=trim(filename), status="replace", action="write")
						else
							open (unit=1, file=trim(filename), status="old", action="write", position="append")
						endif
					endif
				endif

				do idim=1, 1 !ndim
					write (1,*) "zone"
					pstart = 1
					do iphs = 1, nphases
						pend = pstart + phase_array(iphs)%npart- 1
						do m = pstart, pend

!							if (vel_var(iphs,idim)>post_small.and.acc_var(iphs,idim)>post_small) then
								write (1,"(3D15.7)") acc_fluc(m,idim)/acc_var(iphs,idim), vel_fluc(m,idim)/vel_var(iphs,idim), &
									&	acc_fluc(m,idim)*vel_fluc(m,idim) / acc_var(iphs,idim)/vel_var(iphs,idim)
!							endif
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

				filename = trim(run_name)//"_AiVi_avg"
				if (from_post) then
					filename=trim(filename)//"_post.dat"
					open (unit=1, file=trim(filename), status="replace", action="write")
					write (1,*) "zone"
				else
					filename=trim(filename)//".dat"
					if (irestart==0.and.first_pass) then
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
				endif

				tmp1 = sqrt( dot_product (acc_var(1,:), vel_var(1,:)) )
				if (from_post) then
					write (1,"(6f10.3,2D15.7)") dbydx, lybyd, maxvolfrac, re, rhos/rhof, coeff_rest, source/tmp1, dissipation/tmp1
				else
					write (1,"(3D15.7)") t/dia_phys*umeanslip/(1-maxvolfrac), source/tmp1, dissipation/tmp1
				endif
				close (1)
				write (*,"(1A,2d15.7)") "SOURCE, DISSIPATION = ", source/tmp1, dissipation/tmp1
			endif
		else
			if(I_AM_NODE_ZERO)then
				nvar  = 8
				nvar1 = 6
				nvar2 = 2

				filename = "_AiVi_avg.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line)
			endif 
		endif
		endif
	end subroutine AiVi

	subroutine interstitial_dist(part_pos, part_rad, int_dist, npart)
		implicit none
		integer, intent(in) :: npart
		real(prcn), intent(in) :: part_pos(npart,ndim), part_rad(npart)
		real(prcn), intent(out) :: int_dist
		integer :: ibody, jbody, idim, imis, ivar, nvar
		integer, allocatable :: neighb(:)
		real(prcn), allocatable :: length(:)
		real(prcn) :: sep, rad_reg, dist(ndim), distmod, nom, denom, weight
		character*50 filename
		real(prcn), allocatable :: var(:,:), avr_var(:), var_var(:)

		real(prcn), allocatable :: gofr(:), rad_bin(:)
		integer :: nrbins, pair



return

		write (*,*) "COMPUTING INTERSTITIAL DISTANCE ..."

		if (.not.post_no_flow_mem_alloc) then
			allocate(length(npart), neighb(npart))
			length = zero
			neighb = 0

			write (*,"(1A,1f5.2)") "RAD_FACTOR = ", rad_factor
			write (*,"(1A,1f5.2)") "POWER      = ", power

			sep = my/2.0
#if 0
			nom   = zero
			denom = zero
			do ibody=1, npart
				rad_reg = part_rad(ibody)*2*rad_factor
				do jbody=1, npart
					if (jbody/=ibody) then
						dist(:) = abs(part_pos(ibody,:)-part_pos(jbody,:))
						do idim=1, ndim
							if (dist(idim)>sep) dist(idim) = my - dist(idim)
						enddo
				
						distmod = dsqrt(dot_product(dist,dist))

						if (distmod-2*part_rad(ibody)<small_number) then
							write (*,*) "CONTACT DETECTED"
							write (*,*) distmod
!						else
						elseif (distmod<=rad_reg) then
							distmod = distmod - (part_rad(ibody)+part_rad(jbody))
							weight = 1./distmod**power

							nom   = nom   + distmod * weight
							denom = denom + weight
!							nom   = nom   + 1./(distmod-2*part_rad(ibody))
!							denom = denom + 1./(distmod-2*part_rad(ibody))**2
!!!							length(ibody) = length(ibody)+distmod
							neighb(ibody) = neighb(ibody) + 1
						endif
					endif
				enddo
				if (denom>small_number) length(ibody) = nom/denom
!				if (neighb(ibody)>0) length(ibody) = length(ibody)/neighb(ibody)
			enddo
!			length = length * dy/ dia_phys

			int_dist = nom/denom * dia_phys / dbydx

			write (*,"(A,2D15.7)") "NOM, DENOM, PAIRS = ", nom, denom
			write (*,*) "INTERSTITIAL DISTANCE = ", int_dist


!			call screen_separator(30,'^')
!			write (*,*) "WRITING THE SCATTER PLOT"

!			filename = trim(run_name)//"_int_dist_scatter.dat"
!			open(unit=1, file=trim(filename), status="replace")
!			do ibody=1,nbody
!				write (1,"(1F6.2,1I4,1D15.7)") maxvolfrac, neighb(ibody), length(ibody)
!			enddo
!			close (1)

!			int_dist = zero
!			denom = zero
!			do ibody=1, nbody
!				if (length(ibody)>small_number) then
!					int_dist = int_dist + length(ibody)
!					denom = denom + 1.
!				endif
!			enddo
!			int_dist = int_dist/denom

#endif
			nom   = zero
			denom = zero
			pair = 0
			do ibody=1, npart-1
				rad_reg = part_rad(ibody)*2*rad_factor
!				nom   = zero
!				denom = zero
				do jbody=ibody+1, npart
!					if (jbody/=ibody) then
						dist(:) = abs(part_pos(ibody,:)-part_pos(jbody,:))
						do idim=1, ndim
							if (dist(idim)>sep) dist(idim) = my - dist(idim)
						enddo
				
						distmod = dsqrt(dot_product(dist,dist))

						if (distmod - (part_rad(ibody)+part_rad(jbody)) < small_number) then
							write (*,*) "CONTACT DETECTED"
							write (*,*) distmod
!						else
						elseif (distmod<=rad_reg) then
							distmod = distmod - (part_rad(ibody)+part_rad(jbody))
							weight = 1./distmod**power

							nom    = nom + distmod*weight
							denom  = denom + weight
!							neighb(ibody) = neighb(ibody) + 1
						else
							pair = pair + 1
						endif
!					endif
				enddo
!				if (denom>small_number) length(ibody) = nom/denom
!				if (neighb(ibody)>0) length(ibody) = length(ibody)/neighb(ibody)
			enddo
			int_dist = nom/denom
			int_dist = int_dist * dia_phys / dbydx

!			call screen_separator(30,'^')

			write (*,"(A,2D15.7,1I10)") "NOM, DENOM, PAIRS = ", nom, denom, pair
			write (*,*) "INTERSTITIAL DISTANCE = ", int_dist
!			filename = trim(run_name)//"_int_dist.dat"
!			open(unit=1, file=trim(filename), status="replace")
!			write (1,"(1D15.7)") int_dist
!			close (1)
!			call screen_separator(30,'I')

!			nrbins = 200 
!			allocate(gofr(nrbins), rad_bin(nrbins))
			
!			CALL calculate_gofr_homog(nbody,part_pos(1:nbody,1:3), my, mxf, nrbins, .true., gofr(1:nrbins), rad_bin(1:nrbins))
			
		else
			if(I_AM_NODE_ZERO)then
				nvar = 1
				allocate(var(nmis,nvar))
	
				do imis=1, nmis
					if (imis==1) then
						filename="MIS1_int_dist.dat"
					elseif (imis==2) then
						filename="MIS2_int_dist.dat"
					elseif (imis==3) then
						filename="MIS3_int_dist.dat"
					elseif (imis==4) then
						filename="MIS4_int_dist.dat"
					elseif (imis==5) then
						filename="MIS5_int_dist.dat"
					elseif (imis==6) then
						filename="MIS6_int_dist.dat"
					elseif (imis==7) then
						filename="MIS7_int_dist.dat"
					elseif (imis==8) then
						filename="MIS8_int_dist.dat"
					endif
					
					write (*,*) "READING FROM "//filename
					open (unit=1,file=trim(filename),status="old",action="read")
					read (1,"(1D15.7)") var(imis,:)
					close (1)
				enddo

				allocate(avr_var(nvar),var_var(nvar))
				avr_var = 0d0
				var_var = 0d0

				do ivar=1, nvar
					avr_var(ivar) = sum(var(:,ivar))
				enddo
				avr_var(:) = avr_var(:)/nmis

				do ivar=1, nvar
					do imis=1, nmis
						var_var(ivar) = var_var(ivar) + (var(imis,ivar)-avr_var(ivar))**2
					enddo
				enddo
				var_var(:) = var_var(:)/nmis
				var_var(:) = sqrt(var_var(:))


				filename = "total_int_dist.dat"
				open (unit=1,file=trim(filename),status="replace",action="write")
				write (1,"(2D15.7)") avr_var(:), var_var(:)
				close (1)
			endif 
		endif
	end subroutine interstitial_dist


	subroutine int_dist_gofr(nbins, radbin, gofr, int_dist)
		implicit none
		integer :: nbins
		real(prcn) :: radbin(nbins), gofr(nbins), int_dist

		integer :: i
		real(prcn) :: dist, weight, nom, denom, partbin(nbins), dr, rad_reg


		dr = (radbin(2)-radbin(1)) * lybyd
		rad_reg = dia_phys*rad_factor

		nom = zero
		denom = zero
		do i=1, nbins
			if (gofr(i)>0) then
				partbin(i) = gofr(i) * 4. * pi * (radbin(i)*lybyd)**2 * dr * nbody/2 * nbody / lybyd**3
				dist = radbin(i)*lybyd
				if (dist<=rad_reg) then
					dist = dist - dia_phys
					weight = 1. / dist**power
					nom = nom + partbin(i) * dist * weight
					denom = denom + partbin(i) * weight
				endif
			endif
		enddo

		int_dist = nom / denom

		write (*,"(A,2D15.7)") "NOM, DENOM = ", nom, denom
		write (*,"(A,1D15.7)") "INTERSTITIAL DISTANCE = ", int_dist
	end subroutine int_dist_gofr


	subroutine uiui_correlation
		use postproc_funcs
!		use init_turb, only : gener_filename
		implicit none
		integer :: i, j, k, l, i1, i2, j1, j2, k1, k2, ii1, ii2, ijk1, ijk2, imis, idim, nvar, nvar1, nvar2, nline, ivar, ibin
		real(prcn) :: dr, rmax, var_lint, avr_lint, sd_lint, err_bar, confint, rij_mag
		real(prcn), allocatable :: uiui_par(:), uiui_per0(:), uiui_per1(:), uiui_per2(:)
		real(prcn), allocatable :: uiui_par_s(:), uiui_per0_s(:), uiui_per1_s(:), uiui_per2_s(:)
		real(prcn), allocatable :: uiui_par_fs(:), uiui_per0_fs(:), uiui_per1_fs(:), uiui_per2_fs(:)

		integer(8), allocatable :: num_bins(:), num_bins_fs(:), num_bins_s(:)
		real(prcn) :: ui(ndim), uj(ndim), u1i(ndim), u1j(ndim), u2i(ndim), u2j(ndim), u1i_mag, u1j_mag
		real(prcn) :: cpu1, cpu0

		integer :: iproc, id_s, id_r, node_num_s, node_num_r, s_size
		real(prcn), allocatable, dimension(:,:,:,:) :: urecv

		logical, allocatable :: fluid_atijk2(:,:,:)

		character*50 filename1
		character*5 turn
		character*1 rank_string
		integer :: proc_start, ijk_start, count_out, count=0
		integer :: strlen
		logical :: filexist, finish

		rmax  = doml(2) / 2

		allocate(uiui_par(0:nbins), uiui_per0(0:nbins), uiui_per1(0:nbins), uiui_per2(0:nbins))
		allocate(num_bins(0:nbins))

		if (imove==1.or.move_particles) then
			allocate(uiui_par_s(0:nbins), uiui_per0_s(0:nbins), uiui_per1_s(0:nbins), uiui_per2_s(0:nbins))
			allocate(uiui_par_fs(0:nbins), uiui_per0_fs(0:nbins), uiui_per1_fs(0:nbins), uiui_per2_fs(0:nbins))
			allocate(num_bins_fs(0:nbins), num_bins_s(0:nbins))
		endif
		if (allocated(rad_bin)) deallocate(rad_bin)
		allocate(rad_bin(0:nbins))

		dr = rmax / nbins

		id_r = myid
		id_s = myid

		do ibin=0, nbins
			rad_bin(ibin) = ibin * dr
		enddo

		uiui_par  = zero
		uiui_per0 = zero
		uiui_per1 = zero
		uiui_per2 = zero
		num_bins     = zero

		if (imove==1.or.move_particles) then
			uiui_par_s  = zero
			uiui_per0_s = zero
			uiui_per1_s = zero
			uiui_per2_s = zero

			uiui_par_fs  = zero
			uiui_per0_fs = zero
			uiui_per1_fs = zero
			uiui_per2_fs = zero

			num_bins_fs  = zero
			num_bins_s   = zero
		endif

		if (.not.post_no_flow_mem_alloc) then
			call screen_separator(30,'I')
			if (I_AM_NODE_ZERO) write (*,*) "GENERATING THE PARALLEL AND PERPENDICULAR CORRELATIONS"

			call initialize_gridvertindex (nx, my, mz)

			call GENER_FILENAME(filename1,trim(run_name)//"_corr_res.rst")
			inquire(file=filename1, exist=filexist)

			if (filexist) then
				call restart_in
			else
				proc_start = 0
				ijk_start  = 0
			endif

			do iproc=proc_start, nproc/2
				count_out = 0
				first_pass = .true.
				do ijk1 = ijk_start+1, nvert, skip_num
					count_out = count_out+1
					call ind1t3(ijk1,i1,j1,k1)

					if (.not.((imove==1.or.move_particles) .or. fluid_atijk(i1,j1,k1))) goto 20

					if (iproc==0) then
						do ijk2 = ijk1, nvert
							call ind1t3(ijk2,i2,j2,k2)						
!							if (fluid_atijk(i2,j2,k2)) call calc_correlation(ubcp(1:nx,:,:,:), ubcp(1:nx,:,:,:), i1, j1, k1, i2, j2, k2)
							if ((imove==1.or.move_particles) .or. fluid_atijk(i2,j2,k2)) then
								call calc_correlation(ubcp(i1,j1,k1,1:ndim), ubcp(i2,j2,k2,1:ndim), fluid_atijk(i1,j1,k1), fluid_atijk(i2,j2,k2), i1, j1, k1,  i2, j2, k2)
							endif
						enddo
#if PARALLEL
					else
						if (first_pass) then
							id_s = myid+iproc
							id_r = myid-iproc

							if (id_s>nproc-1) id_s = id_s-nproc
							if (id_r<0)       id_r = id_r+nproc

							node_num_s = nx*my*mz
							node_num_r = (ends(id_r)-starts(id_r)+1)*my*mz

							if (allocated(urecv)) deallocate(urecv,fluid_atijk2)
							allocate(urecv(starts(id_r):ends(id_r), my, mz, ndim))
							allocate(fluid_atijk2(starts(id_r):ends(id_r), my, mz))

							CALL MPI_SENDRECV(ubcp(1:nx,:,:,:), node_num_s*ndim, MPI_DOUBLE_PRECISION, id_s, myid, urecv, node_num_r*ndim, MPI_DOUBLE_PRECISION, id_r, id_r, decomp_group, status, err_code)

							CALL MPI_SENDRECV(fluid_atijk(1:nx,:,:), node_num_s, MPI_LOGICAL, id_s, id_s, fluid_atijk2, node_num_r, MPI_LOGICAL, id_r, myid, decomp_group, status, err_code)

							first_pass = .false.
						endif

						if (iproc==nproc/2 .and. myid+1>nproc/2) goto 10

						do k2=1, mz
							do j2=1, my
								do i2=starts(id_r), ends(id_r)
!									if (fluid_atijk2(i2,j2,k2)) call calc_correlation(ubcp(1:nx,:,:,:), urecv, starts(myid)+i1-1, j1, k1, i2, j2, k2)
									if ((imove==1.or.move_particles) .or. fluid_atijk(i2,j2,k2)) then
										call calc_correlation(ubcp(i1,j1,k1,1:ndim), urecv(i2,j2,k2,1:ndim), fluid_atijk(i1,j1,k1), fluid_atijk2(i2,j2,k2), starts(myid)+i1-1, j1, k1,  i2, j2, k2)
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
						write (*,"(4(1a,1i))") " NODE ", ijk1, " OUT OF ", nx*my*mz, " BETWEEN PROCs ", myid, " AND ", id_r
#else
						write (*,"(2(1a,1i))") " NODE ", ijk1, " OUT OF ", nx*my*mz
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
#if PARALLEL
					BROADCAST_DOUBLE(cpu1,1,NODE_ZERO,decomp_group)
#endif
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

			deallocate (uiui_par, uiui_per0, uiui_per1, uiui_per2)
			deallocate (rad_bin, num_bins)

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
				call mis_average(nvar, nvar1, nvar2, filename1, line)
			endif 
		endif

		if (I_AM_NODE_ZERO) call screen_separator(30,'I')

	contains

		subroutine calc_correlation(u1, u2, fluid1, fluid2, i1, j1, k1, i2, j2, k2)
			implicit none
			integer, intent(in) :: i1, i2, j1, j2, k1, k2
			real(prcn), intent(in) :: u1(ndim), u2(ndim)
			logical, intent(in) :: fluid1, fluid2
			real(prcn) :: r, rij(ndim), tmp(ndim+1), vec1(ndim), vec2(ndim)
!			integer :: imax1, imax2, jmax1, jmax2, kmax1, kmax2
			integer :: body1, body2

!			imax1 = size(u1,1)
!			jmax1 = size(u1,2)
!			kmax1 = size(u1,3)

!			imax2 = size(u2,1)
!			jmax2 = size(u2,2)
!			kmax2 = size(u2,3)

			rij(1) = abs(i2-i1) * dx
			rij(2) = abs(j2-j1) * dy
			rij(3) = abs(k2-k1) * dz

			if (.not.fluid1) call find_body(i1,j1,k1,body1, rmax)
			if (.not.fluid2) call find_body(i2,j2,k2,body2, rmax)

!			if (.not.fluid1) write (*,*) "body : ", body1
!			if (.not.fluid2) write (*,*) "body : ", body2

			do idim=1, ndim
				if (rij(idim) > rmax) rij(idim) = doml(idim)-rij(idim)
			enddo

			r = sqrt(rij(1)**2 + rij(2)**2 + rij(3)**2) 
			if (r<=rmax) then
				ibin = int(r/rmax * nbins) +1
				if (ibin>nbins) ibin = nbins
				if (i1==i2.and.j1==j2.and.k1==k2) ibin = 0

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

				do idim=1, ndim
					tmp(idim) = vec1(idim) * vec2(idim)
				enddo
				tmp(ndim+1) = tmp(2)+tmp(3)

				if (fluid1.and.fluid2) then
					num_bins(ibin)  = num_bins(ibin)  + 1
					uiui_par(ibin)  = uiui_par(ibin)  + tmp(1)
					uiui_per1(ibin) = uiui_per1(ibin) + tmp(2)
					uiui_per2(ibin) = uiui_per2(ibin) + tmp(3)
					uiui_per0(ibin) = uiui_per0(ibin) + tmp(4)
				elseif (fluid1.or.fluid2) then
					num_bins_fs(ibin)  = num_bins_fs(ibin)  + 1
					uiui_par_fs(ibin)  = uiui_par_fs(ibin)  + tmp(1)
					uiui_per1_fs(ibin) = uiui_per1_fs(ibin) + tmp(2)
					uiui_per2_fs(ibin) = uiui_per2_fs(ibin) + tmp(3)
					uiui_per0_fs(ibin) = uiui_per0_fs(ibin) + tmp(4)
				else
					num_bins_s(ibin)  = num_bins_s(ibin)  + 1
					uiui_par_s(ibin)  = uiui_par_s(ibin)  + tmp(1)
					uiui_per1_s(ibin) = uiui_per1_s(ibin) + tmp(2)
					uiui_per2_s(ibin) = uiui_per2_s(ibin) + tmp(3)
					uiui_per0_s(ibin) = uiui_per0_s(ibin) + tmp(4)
				endif
			endif
		end subroutine calc_correlation

		subroutine output
			implicit none
			character*50 filename
			real(prcn), dimension(0:nbins) :: u1out, u2out, u3out, u4out
			real(prcn), dimension(0:nbins) :: u1out_s, u2out_s, u3out_s, u4out_s
			real(prcn), dimension(0:nbins) :: u1out_fs, u2out_fs, u3out_fs, u4out_fs
			integer(8), dimension(0:nbins) :: num_binsout, num_binsout_fs, num_binsout_s
			integer :: i, j

#if PARALLEL
			GLOBAL_DOUBLE_SUM(uiui_par, u1out, nbins+1, decomp_group)
			GLOBAL_DOUBLE_SUM(uiui_per1, u2out, nbins+1, decomp_group)
			GLOBAL_DOUBLE_SUM(uiui_per2, u3out, nbins+1, decomp_group)
			GLOBAL_DOUBLE_SUM(uiui_per0, u4out, nbins+1, decomp_group)
			call mpi_allreduce(num_bins, num_binsout, nbins+1, mpi_integer8, mpi_sum, decomp_group, err_code)

			if (imove==1.or.move_particles) then
				GLOBAL_DOUBLE_SUM(uiui_par_s, u1out_s, nbins+1, decomp_group)
				GLOBAL_DOUBLE_SUM(uiui_per1_s, u2out_s, nbins+1, decomp_group)
				GLOBAL_DOUBLE_SUM(uiui_per2_s, u3out_s, nbins+1, decomp_group)
				GLOBAL_DOUBLE_SUM(uiui_per0_s, u4out_s, nbins+1, decomp_group)
				call mpi_allreduce(num_bins_s, num_binsout_s, nbins+1, mpi_integer8, mpi_sum, decomp_group, err_code)

				GLOBAL_DOUBLE_SUM(uiui_par_fs, u1out_fs, nbins+1, decomp_group)
				GLOBAL_DOUBLE_SUM(uiui_per1_fs, u2out_fs, nbins+1, decomp_group)
				GLOBAL_DOUBLE_SUM(uiui_per2_fs, u3out_fs, nbins+1, decomp_group)
				GLOBAL_DOUBLE_SUM(uiui_per0_fs, u4out_fs, nbins+1, decomp_group)
				call mpi_allreduce(num_bins_fs, num_binsout_fs, nbins+1, mpi_integer8, mpi_sum, decomp_group, err_code)
			endif
#else
			u1out = uiui_par
			u2out = uiui_per1
			u3out = uiui_per2
			u4out = uiui_per0
			num_binsout = num_bins

			if (imove==1.or.move_particles) then
				u1out_s = uiui_par_s
				u2out_s = uiui_per1_s
				u3out_s = uiui_per2_s
				u4out_s = uiui_per0_s
				num_binsout_s = num_bins_s

				u1out_fs = uiui_par_fs
				u2out_fs = uiui_per1_fs
				u3out_fs = uiui_per2_fs
				u4out_fs = uiui_per0_fs
				num_binsout_fs = num_bins_fs
			endif
#endif
			if (I_AM_NODE_ZERO) then
				do ibin=0, nbins
					if (num_binsout(ibin)>0) then
						u1out(ibin) = u1out(ibin) / num_binsout(ibin)
						u2out(ibin) = u2out(ibin) / num_binsout(ibin)
						u3out(ibin) = u3out(ibin) / num_binsout(ibin)
						u4out(ibin) = u4out(ibin) / num_binsout(ibin)
					endif

					if (imove==1.or.move_particles) then
						if (num_binsout_s(ibin)>0) then
							u1out_s(ibin) = u1out_s(ibin) / num_binsout_s(ibin)
							u2out_s(ibin) = u2out_s(ibin) / num_binsout_s(ibin)
							u3out_s(ibin) = u3out_s(ibin) / num_binsout_s(ibin)
							u4out_s(ibin) = u4out_s(ibin) / num_binsout_s(ibin)
						endif

						if (num_binsout_fs(ibin)>0) then
							u1out_fs(ibin) = u1out_fs(ibin) / num_binsout_fs(ibin)
							u2out_fs(ibin) = u2out_fs(ibin) / num_binsout_fs(ibin)
							u3out_fs(ibin) = u3out_fs(ibin) / num_binsout_fs(ibin)
							u4out_fs(ibin) = u4out_fs(ibin) / num_binsout_fs(ibin)
						endif
					endif
				enddo

				do ibin=nbins ,0, -1
					if (u1out(0)>small_number) u1out(ibin) = u1out(ibin) / u1out(0)
					if (u2out(0)>small_number) u2out(ibin) = u2out(ibin) / u2out(0)
					if (u3out(0)>small_number) u3out(ibin) = u3out(ibin) / u3out(0)
					if (u4out(0)>small_number) u4out(ibin) = u4out(ibin) / u4out(0)

					if (imove==1.or.move_particles) then
						if (u1out_s(0)>small_number) u1out_s(ibin) = u1out_s(ibin) / u1out_s(0)
						if (u2out_s(0)>small_number) u2out_s(ibin) = u2out_s(ibin) / u2out_s(0)
						if (u3out_s(0)>small_number) u3out_s(ibin) = u3out_s(ibin) / u3out_s(0)
						if (u4out_s(0)>small_number) u4out_s(ibin) = u4out_s(ibin) / u4out_s(0)

						if (u1out_fs(1)>small_number) u1out_fs(ibin) = u1out_fs(ibin) / u1out_fs(1)
						if (u2out_fs(1)>small_number) u2out_fs(ibin) = u2out_fs(ibin) / u2out_fs(1)
						if (u3out_fs(1)>small_number) u3out_fs(ibin) = u3out_fs(ibin) / u3out_fs(1)
						if (u4out_fs(1)>small_number) u4out_fs(ibin) = u4out_fs(ibin) / u4out_fs(1)
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

				do ibin=0,nbins
!					if(num_binsout(ibin)>0) then
						if (imove==1.or.move_particles) then
							write (1,'(1f10.4,12d15.7)') rad_bin(ibin), u1out(ibin), u2out(ibin), u3out(ibin), u4out(ibin), u1out_s(ibin), u2out_s(ibin), u3out_s(ibin), u4out_s(ibin), u1out_fs(ibin), u2out_fs(ibin), u3out_fs(ibin), u4out_fs(ibin)
						else
							write (1,'(1f10.4,4d15.7)') rad_bin(ibin), u1out(ibin), u2out(ibin), u3out(ibin), u4out(ibin)
						endif
!					endif
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
					if (i/=node_zero) SEND_STRING(filename2, strlen, i, 0, 1, decomp_group)
				enddo
			else
				RECV_STRING(filename2, strlen, node_zero, 0, 1, decomp_group, status)
			endif
#else
			filename2 = trim(run_name)//"_corr_res"//trim(turn)//".rst"
#endif	

			if (I_AM_NODE_ZERO) write (*,*) "RESTARTING GENERATION OF AUTOCORRELATION FUNCTION"

			open (unit=9996,file=filename2,status="old",action="read",form="unformatted")

			if (imove==1.or.move_particles) then
				read (9996) proc_start, ijk_start, uiui_par, uiui_per1, uiui_per2, uiui_per0, num_bins, uiui_par_s, uiui_per1_s, uiui_per2_s, uiui_per0_s, num_bins_s, uiui_par_fs, uiui_per1_fs, uiui_per2_fs, uiui_per0_fs, num_bins_fs
			else
				read (9996) proc_start, ijk_start, uiui_par, uiui_per1, uiui_per2, uiui_per0, num_bins
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
					if (i/=node_zero) SEND_STRING(filename, strlen, i, 0, 1, decomp_group)
				enddo
			else
				RECV_STRING(filename, strlen, node_zero, 0, 1, decomp_group, status)
			endif
#else
			filename = trim(run_name)//"_corr_res"//trim(turn)//".rst"
#endif
			open (unit=9998,file=trim(filename),status="replace",action="write",form="unformatted")
			if (imove==1.or.move_particles) then
				write (9998) iproc, ijk1, uiui_par, uiui_per1, uiui_per2, uiui_per0, num_bins, uiui_par_s, uiui_per1_s, uiui_per2_s, uiui_per0_s, num_bins_s, uiui_par_fs, uiui_per1_fs, uiui_per2_fs, uiui_per0_fs, num_bins_fs
			else
				write (9998) iproc, ijk1, uiui_par, uiui_per1, uiui_per2, uiui_per0, num_bins
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

#if 0
	subroutine compute_gofr
		use postproc_funcs
		implicit none
		integer :: nbin, gofunit, gofavgunit, imis, j
		real(prcn), allocatable :: gofr(:), radbin(:), gofravg(:)
		real(prcn) :: tmp
		character*50 filename
		logical :: rescaling

		nbin = 100
		rescaling = .true.
		allocate(gofr(nbin),radbin(nbin),gofravg(nbin))

		gofunit = getnewunit(minunitno,maxunitno)



		radbin = zero
		gofr = zero
		gofravg = zero



		if (.not.post_no_flow_mem_alloc) then
			CALL calculate_gofr_homog(nbody,xc(1:nbody,1:3), my, mxf, nbin, rescaling, gofr, radbin)

			OPEN(gofunit, file=TRIM(RUN_NAME)//'_gofr.dat', form='formatted', status="replace")
			do j=1,nbin
				write (gofunit,'(3D15.7)') radbin(j), radbin(j)*lybyd, gofr(j)
			end do
			close (gofunit, status = "keep")
		else
			gofravg = zero

			if(I_AM_NODE_ZERO)then
				do imis=1, nmis
					if (imis==1) then
						filename="MIS1_gofr.dat"
					elseif (imis==2) then
						filename="MIS2_gofr.dat"
					elseif (imis==3) then
						filename="MIS3_gofr.dat"
					elseif (imis==4) then
						filename="MIS4_gofr.dat"
					elseif (imis==5) then
						filename="MIS5_gofr.dat"
					elseif (imis==6) then
						filename="MIS6_gofr.dat"
					elseif (imis==7) then
						filename="MIS7_gofr.dat"
					elseif (imis==8) then
						filename="MIS8_gofr.dat"
					endif

					write (*,*) "READING FROM "//filename
					open (unit=gofunit,file=trim(filename),status="old",action="read")
					do j=1,nbin
						read (gofunit,'(3D15.7)') radbin(j), tmp, gofr(j)
						gofravg(j) = gofravg(j) + gofr(j)
					enddo
					close (gofunit)
				enddo

				gofravg = gofravg/nmis

				gofavgunit = getnewunit(minunitno,maxunitno)
				open(gofavgunit, file='gofr_avg.dat', form='formatted', status="replace")
				do j=1,nbin
					write (gofavgunit,'(3D15.7)') radbin(j), radbin(j)*lybyd, gofravg(j)
				end do
				close (gofavgunit, status = "keep")
			endif
		endif
	end subroutine compute_gofr
#endif

#if 0
	subroutine slice_umean_comp
		implicit none
		integer :: iconfig, nconfig, islice
		integer :: i,j,k
		character*50 filename

!		call calc_pgrad

		nconfig = 0
		do i=1, my/2
			if (mod(my,i)==0) then
				nconfig = nconfig + 1
			endif
		enddo

		call screen_separator(80,'I')
		write (*,*) "NUMBER OF CONFIGURATIONS = ", nconfig

		allocate(slice(nconfig))

		iconfig = 0
		do i=1, my/2
			if (mod(my,i)==0) then
				iconfig = iconfig + 1
				slice(iconfig)%num = my/i
				slice(iconfig)%dy  = i

				allocate(slice(iconfig)%ind(slice(iconfig)%num))
				allocate(slice(iconfig)%fluidnum(slice(iconfig)%num))
				allocate(slice(iconfig)%solidnum(slice(iconfig)%num))
				allocate(slice(iconfig)%ufmean(slice(iconfig)%num,ndim))
				allocate(slice(iconfig)%mpg(slice(iconfig)%num,ndim))
				allocate(slice(iconfig)%areaf(slice(iconfig)%num))

				slice(iconfig)%ind = 0
				slice(iconfig)%fluidnum = 0
				slice(iconfig)%solidnum = 0
				slice(iconfig)%ufmean   = 0d0
				slice(iconfig)%mpg      = zero
				slice(iconfig)%areaf    = zero

				do j=1, slice(iconfig)%num
					slice(iconfig)%ind(j) = (j-1)*slice(iconfig)%dy+1
				enddo

				call screen_separator(10,'^')
				write (*,*) "CONFIG #    = ", iconfig
				write (*,*) "DY          = ", slice(iconfig)%dy
				write (*,*) "# OF SLICES = ", slice(iconfig)%num
				call screen_separator(10,'-')
			endif
		enddo

		do k=1, my
			do j=1, my
				do i=1, nx
!					write (*,*) i,j,k
					do iconfig=1, nconfig
						call find_slice_ind(j,iconfig,islice)

						if (fluid_atijk(i,j,k)) then
							slice(iconfig)%fluidnum(islice) = slice(iconfig)%fluidnum(islice) + 1
							slice(iconfig)%ufmean(islice,:) = slice(iconfig)%ufmean(islice,:) + ubcp(i,j,k,:)
							slice(iconfig)%mpg(islice,:) = slice(iconfig)%mpg(islice,:) + ppr(i,j,k,:)
							slice(iconfig)%areaf(islice) = slice(iconfig)%areaf(islice) + 1
						else
							slice(iconfig)%solidnum(islice) = slice(iconfig)%solidnum(islice) + 1
						endif
					enddo
				enddo
			enddo
		enddo

!do iconfig=1, nconfig
!do islice=1, slice(iconfig)%num
!write (*,*) slice(iconfig)%ufmean(islice,:)

!read (*,*) 


		do iconfig=1, nconfig
			do islice=1, slice(iconfig)%num
				slice(iconfig)%ufmean(islice,:) = slice(iconfig)%ufmean(islice,:) / slice(iconfig)%fluidnum(islice) / umeanslip
				slice(iconfig)%mpg(islice,:) = -slice(iconfig)%mpg(islice,:) / slice(iconfig)%fluidnum(islice) / (umeanslip**2/2/dia_phys)
				slice(iconfig)%areaf(islice) = slice(iconfig)%areaf(islice) / my**2
			enddo
		enddo



		filename = trim(run_name)//"_UFMEANY.dat"
		open (unit=1, file=trim(filename), status="replace", action="write")

		write (1,*) "variables=y,u1,u2,u3,p1,p2,p3,a"
		do iconfig=1, nconfig
			write (1,*) 'zone t= "', slice(iconfig)%dy, '"'
			do islice=1, slice(iconfig)%num
				write (1,"(1I6,7D15.7)") (islice-1)*slice(iconfig)%dy+slice(iconfig)%dy/2, slice(iconfig)%ufmean(islice,:), slice(iconfig)%mpg(islice,:), &
							& 	slice(iconfig)%areaf(islice)
			enddo
		enddo
		close (1)
	end subroutine slice_umean_comp


	subroutine find_slice_ind(j,iconfig,ind)
		implicit none
		integer, intent(in)  :: j, iconfig
		integer, intent(out) :: ind

		integer :: i
	
		do i=1, slice(iconfig)%num
			if (j>=slice(iconfig)%ind(i).and.j<slice(iconfig)%ind(i)+slice(iconfig)%dy) then
				ind = i
!				write (*,*) iconfig, i
				return
			endif
		enddo

		if (i==slice(iconfig)%num+1) then
			write (*,*) "CANNOT	FIND THE SLICE"
			write (*,*) iconfig
			stop
		endif
	end subroutine find_slice_ind
#endif


	subroutine compute_sijsij
		implicit none
	
		real(prcn), allocatable :: dissip(:,:,:), fluc(:,:,:)
		real(prcn), allocatable, dimension(:) :: pres_part, visc_part, int_tke, pres_part_fluc, visc_part_fluc, int_tke_fluc

		real(prcn), allocatable :: total_force(:,:), int_mom(:,:), mix_total_force(:), mix_int_mom(:)
		real(prcn), allocatable :: total_force_mag(:), int_mom_mag(:)
		real(prcn) :: mix_total_force_mag, mix_int_mom_mag
 
		real(prcn) :: mix_pres_part, mix_visc_part, mix_int_tke, mix_pres_part_fluc, mix_visc_part_fluc, mix_int_tke_fluc
		real(prcn) :: sijsij_loc, sijsij_far, sijsij_p, sijsij_far_loc, sijsij_p_loc
		character*3 mis_str
		character*50 filename
	
		real(prcn) :: tmp1, tmp2
		integer :: i, j, k, m, idim, iphs, part_start, part_end
		integer :: im, ip, jm, jp, km, kp, ii, jj, kk, iii, jjj, kkk
		integer :: dim1, dim2, dim3
		logical :: neighb_insolid

		integer :: nvar, nvar1, nvar2
		character*5 mis_string

		write (*,*) "IN COMPUTE_SIJSIJ"
		if (.not.post_no_flow_mem_alloc) then

			call compute_mix_mean_vel

!		allocate (dissip(nx,my,mz),fluc(nx,my,mz))
!		dissip = 0d0
!		fluc   = 0d0

			if (near_particle_checked==.false.) then
				call near_particle_region
				call find_fluid_plus
				near_particle_checked = .true.
			endif

			sijsij_loc     = 0d0
			sijsij_p_loc   = 0d0
			sijsij_far_loc = 0d0
			! generating thetaij
			if (I_AM_NODE_ZERO) write (*,*) "GENERATING SijSij"
			do i=1, nx
				if (I_AM_NODE_ZERO.and.mod(i,20)==0) write (*,*) " PLANE # = ", i
				do dim1=1, ndim
					do dim2=1, ndim
						call derivative_1st(i,dim1,dim2,ur1)
						call derivative_1st(i,dim2,dim1,ur2)

						do k=1, mz
							do j=1, my
								if (fluid_atijk(i,j,k)) then
									tmp1 = 2*vis*((ur1(j,k)+ur2(j,k))*half)**2
									sijsij_loc = sijsij_loc + tmp1
		
	!								dissip(i,j,k) = dissip(i,j,k) + tmp1
	!								fluc(i,j,k) = dot_product(ubcp(i,j,k,:)-ufmean(:),ubcp(i,j,k,:)-ufmean(:))/2

!									if (fluid_atijkp(i,j,k)) then
!										sijsij_p_loc = sijsij_p_loc + tmp1
!									endif

									if (far_field(i,j,k)) then
										sijsij_far_loc = sijsij_far_loc + tmp1
									endif
								else
	!								dissip(i,j,k) = -(vis*(umeanslip/dia_phys)**2)
	!								fluc(i,j,k)   = -(dot_product(ufmean(:),ufmean(:))/2)
								endif
							enddo
						enddo
					enddo
				enddo
			enddo
	!		dissip = dissip / (vis*(umeanslip/dia_phys)**2)
	!		fluc   = fluc   / (dot_product(ufmean(:)-usmean(:),ufmean(:)-usmean(:))/2)

#if PARALLEL
			sijsij = 0d0
			GLOBAL_DOUBLE_SUM(sijsij_loc,sijsij,1,decomp_group)

!			sijsij_p = 0d0
!			GLOBAL_DOUBLE_SUM(sijsij_p_loc,sijsij_p,1,decomp_group)

			sijsij_far = 0d0
			GLOBAL_DOUBLE_SUM(sijsij_far_loc,sijsij_far,1,decomp_group)
#else
			sijsij     = sijsij_loc
!			sijsij_p   = sijsij_p_loc
			sijsij_far = sijsij_far_loc
#endif

			sijsij_pure= sijsij     / count_fluid
			sijsij     = sijsij     / count_fluid / (vis*(mix_mean_slip_mag/dia_phys)**2)
!			sijsij_p   = sijsij_p   / count_fluid / (vis*(mix_mean_slip_mag/dia_phys)**2)
			sijsij_far = sijsij_far / count_fluid / (vis*(mix_mean_slip_mag/dia_phys)**2)

			allocate(pres_part(nphases), visc_part(nphases), int_tke(nphases))
			allocate(pres_part_fluc(nphases), visc_part_fluc(nphases), int_tke_fluc(nphases))

			allocate(int_mom(nphases,ndim), total_force(nphases,ndim))
			allocate(mix_int_mom(ndim), mix_total_force(ndim))
			allocate(int_mom_mag(nphases), total_force_mag(nphases))

			int_tke     = zero
			pres_part   = zero
			visc_part   = zero

			int_mom     = zero
			total_force = zero

			mix_int_tke   = zero
			mix_pres_part = zero
			mix_visc_part = zero

			do iphs=1, nphases
				part_start = phase_array(iphs)%pstart
				part_end   = phase_array(iphs)%pend

				write (*,*) "IPHASE = " , iphs
				write (*,*) "START, END = " , part_start, part_end

				do m=part_start, part_end
					do idim=1, ndim
						tmp1 = ufmean(idim)-phase_array(iphs)%mean_spec_vel(idim)
						pres_part(iphs)   = pres_part(iphs) + tmp1 * pres(m,idim)
						visc_part(iphs)   = visc_part(iphs) + tmp1 * visc(m,idim)
						int_tke(iphs)     = int_tke(iphs)   + tmp1 * (visc(m,idim)+pres(m,idim))

						mix_pres_part   = mix_pres_part + tmp1 * pres(m,idim) !* phase_array(iphs)%volfrac
						mix_visc_part   = mix_visc_part + tmp1 * visc(m,idim) !* phase_array(iphs)%volfrac
						mix_int_tke     = mix_int_tke   + tmp1 * (visc(m,idim)+pres(m,idim)) !* phase_array(iphs)%volfrac

						int_mom(iphs,idim)     = int_mom(iphs,idim)     + (visc(m,idim)+pres(m,idim))
						total_force(iphs,idim) = total_force(iphs,idim) + (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)
					enddo
				enddo
			enddo


usmean1_mag = sqrt(dot_product(phase_array(1)%mean_spec_vel(:),phase_array(1)%mean_spec_vel(:)))
usmean2_mag = sqrt(dot_product(phase_array(2)%mean_spec_vel(:),phase_array(2)%mean_spec_vel(:)))

ufmean_mag = sqrt(dot_product(ufmean(:),ufmean(:)))

mix_mean_vel_mag = (usmean1_mag*phase_array(1)%volfracg+usmean2_mag*phase_array(2)%volfracg)/maxvolfrac

mix_mean_slip_mag = ufmean_mag-mix_mean_vel_mag

write (*,"(4d15.7)") ufmean_mag/mix_mean_slip_mag, mix_mean_vel_mag/mix_mean_slip_mag, usmean1_mag/mix_mean_slip_mag, usmean2_mag/mix_mean_slip_mag

write (*,"(3d15.7)") sqrt(dot_product(mesh_vel,mesh_vel))


			int_tke_fluc   = zero
			pres_part_fluc = zero
			visc_part_fluc = zero

			mix_int_tke_fluc   = zero
			mix_pres_part_fluc = zero
			mix_visc_part_fluc = zero

			if (imove==1) then
				do iphs=1, nphases
					part_start = phase_array(iphs)%pstart
					part_end = phase_array(iphs)%pend

					do m=part_start, part_end
						do idim=1, ndim
							tmp1 = velbdy(m,idim) - phase_array(iphs)%mean_spec_vel(idim)

							pres_part_fluc(iphs) = pres_part_fluc(iphs) - pres(m,idim) * tmp1
							visc_part_fluc(iphs) = visc_part_fluc(iphs) - visc(m,idim) * tmp1
							int_tke_fluc(iphs)   = int_tke_fluc(iphs) - (visc(m,idim)+pres(m,idim)) * tmp1

							mix_pres_part_fluc = mix_pres_part_fluc - pres(m,idim) * tmp1
							mix_visc_part_fluc = mix_visc_part_fluc - visc(m,idim) * tmp1
							mix_int_tke_fluc   = mix_int_tke_fluc   - (visc(m,idim)+pres(m,idim)) * tmp1
						enddo
					enddo
				enddo
			endif

			!^^^ NORMALIZATION ^^^^^^
			do iphs=1, nphases
				int_mom(iphs,:)     = int_mom(iphs,:)     /phase_array(iphs)%npart
				total_force(iphs,:) = total_force(iphs,:) /phase_array(iphs)%npart
			enddo

			do iphs=1, nphases
				int_mom(iphs,:)     = int_mom(iphs,:)     /(3.d0*pi*vis*(mix_mean_slip_mag)*phase_array(iphs)%dia)
				total_force(iphs,:) = total_force(iphs,:) /(3.d0*pi*vis*(mix_mean_slip_mag)*phase_array(iphs)%dia)
			enddo


			mix_int_mom     = zero
			mix_total_force = zero
			do iphs=1, nphases
				mix_int_mom(:)     = mix_int_mom(:)     + int_mom(iphs,:)     * phase_array(iphs)%volfracg/maxvolfrac
				mix_total_force(:) = mix_total_force(:) + total_force(iphs,:) * phase_array(iphs)%volfracg/maxvolfrac
			enddo

			do iphs=1, nphases
				int_mom_mag(iphs)     = sqrt(dot_product(int_mom(iphs,:),int_mom(iphs,:)))
				total_force_mag(iphs) = sqrt(dot_product(total_force(iphs,:),total_force(iphs,:)))
			enddo

			mix_int_mom_mag     = sqrt(dot_product(mix_int_mom(:), mix_int_mom(:)))
			mix_total_force_mag = sqrt(dot_product(mix_total_force(:), mix_total_force(:)))


			do i=1, nphases
				pres_part(i) = pres_part(i) / (vis*(mix_mean_slip_mag/dia_phys)**2)
				visc_part(i) = visc_part(i) / (vis*(mix_mean_slip_mag/dia_phys)**2)
				int_tke(i)   = int_tke(i)   / (vis*(mix_mean_slip_mag/dia_phys)**2)

				pres_part_fluc(i) = pres_part_fluc(i) / (vis*(mix_mean_slip_mag/dia_phys)**2)
				visc_part_fluc(i) = visc_part_fluc(i) / (vis*(mix_mean_slip_mag/dia_phys)**2)
				int_tke_fluc(i)   = int_tke_fluc(i)   / (vis*(mix_mean_slip_mag/dia_phys)**2)
			enddo
			mix_pres_part = mix_pres_part / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) &
										&	/ (vis*(mix_mean_slip_mag/dia_phys)**2)
			mix_visc_part = mix_visc_part / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) &
										&	/ (vis*(mix_mean_slip_mag/dia_phys)**2)
			mix_int_tke   = mix_int_tke   / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) &
										&	/ (vis*(mix_mean_slip_mag/dia_phys)**2)

			mix_pres_part_fluc = mix_pres_part_fluc / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) &
										&	/ (vis*(mix_mean_slip_mag/dia_phys)**2)
			mix_visc_part_fluc = mix_visc_part_fluc / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) &
										&	/ (vis*(mix_mean_slip_mag/dia_phys)**2)
			mix_int_tke_fluc   = mix_int_tke_fluc   / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) &
										&	/ (vis*(mix_mean_slip_mag/dia_phys)**2)

!			do i=1, nphases
!				pres_part_fluc(i) = pres_part_fluc(i) / pres_part(i)
!				visc_part_fluc(i) = visc_part_fluc(i) / visc_part(i)
!				int_tke_fluc(i)   = int_tke_fluc(i)   / int_tke(i)
!			enddo
!			mix_pres_part_fluc = mix_pres_part_fluc / mix_pres_part
!			mix_visc_part_fluc = mix_visc_part_fluc / mix_visc_part
!			mix_int_tke_fluc   = mix_int_tke_fluc   / mix_int_tke

			if (I_AM_NODE_ZERO) then
				filename = trim(run_name)//"_sijsij.dat"
				open (unit=1, file=trim(filename), status="replace")

				if (.not.imove==1) then
					if (nphases==1) then
						write (1,"(4f10.4,9D15.7)") dbydx, lybyd, maxvolfrac, re, &
						& sijsij, sijsij_far, sijsij_far/sijsij, &

						& (int_tke(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke, &
						& (int_tke_fluc(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke_fluc(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke_fluc
					else
						write (1,"(8f10.4,13D15.7)") dbydx, lybyd, maxvolfrac, &
						& (phase_array(iphs)%volfracg, iphs=1,nphases), &
						& (yalpha(iphs),iphs=1,nphases), re,			&
						& sijsij, sijsij_far, sijsij_far/sijsij,		&
						& (int_tke(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke, &
						& (int_tke_fluc(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke_fluc(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke_fluc
					endif
				else
					if (nphases==1) then
						write (1,"(6f10.4,9D15.7)") dbydx, lybyd, maxvolfrac, &
						& re, rhos, coeff_rest,	&
						& sijsij, sijsij_far, sijsij_far/sijsij,	&
						& (int_tke(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke, &
						& (int_tke_fluc(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke_fluc(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke_fluc
					else
						write (1,"(10f10.4,13D15.7)") dbydx, lybyd, maxvolfrac, &
						& (phase_array(iphs)%volfracg, iphs=1,nphases), &
						& (yalpha(iphs), iphs=1,nphases), re,			&
						& rhos, coeff_rest,								&
						& sijsij, sijsij_far, sijsij_far/sijsij,	&
						& (int_tke(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke, &
						& (int_tke_fluc(iphs) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)), iphs=1,nphases),&
						& (int_tke_fluc(iphs) / phase_array(iphs)%npart, iphs=1,nphases),&
						& mix_int_tke_fluc
					endif
				endif
				close (1)

				if (.not.imove==1) then
					filename = trim(run_name)//"_force_comp.dat"
					open (unit=1, file=trim(filename), status="replace")
					if (nphases==1) then
						write (1,"(4f8.4,2D15.7)") dbydx, lybyd, maxvolfrac, re,	&
						& (int_mom_mag(iphs), iphs=1, nphases), (total_force_mag(iphs), iphs=1, nphases)
					else
						write (1,"(8f8.4,6D15.7)") dbydx, lybyd, maxvolfrac, &
						& (phase_array(iphs)%volfracg, iphs=1,nphases),	&
						& (yalpha(iphs), iphs=1,nphases), re,											&
						& (int_mom_mag(iphs), iphs=1, nphases), (total_force_mag(iphs), iphs=1, nphases),&
						& mix_int_mom_mag, mix_total_force_mag
					endif
					close(1)
				endif
			endif
		else
			if(I_AM_NODE_ZERO)then
				if (.not.imove==1) then
					if (nphases==1) then
						nvar  = 13
						nvar1 = 4
						nvar2 = 9
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
				filename = "_sijsij.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line)

				if (.not.imove==1) then
					if (nphases==1) then
						nvar  = 6
						nvar1 = 4
						nvar2 = 2
					else
						nvar  = 14
						nvar1 = 8
						nvar2 = 6
					endif
					filename = "_force_comp.dat"
					call mis_average(nvar, nvar1, nvar2, filename, line)
				endif
			endif 
		endif	
	end subroutine compute_sijsij



#if 0
	subroutine compute_sijsij
		implicit none
	
		real(prcn), allocatable :: dissip(:,:,:), fluc(:,:,:)
		real(prcn) :: pres_part, visc_part, int_tke, pres_part_fluc, visc_part_fluc, int_tke_fluc

		real(prcn), dimension(ndim) :: total_force, int_mom
		real(prcn) :: total_force_mag, int_mom_mag
 
		real(prcn) :: sijsij, sijsij_loc, sijsij_far, sijsij_far_loc
		character*3 mis_str
		character*50 filename
	
		real(prcn) :: tmp1, tmp2
		integer :: i, j, k, m, idim, iphs, part_start, part_end
		integer :: im, ip, jm, jp, km, kp, ii, jj, kk, iii, jjj, kkk
		integer :: dim1, dim2, dim3
		logical :: neighb_insolid

		integer :: nvar, nvar1, nvar2
		character*5 mis_string

		write (*,*) "IN COMPUTE_SIJSIJ"
		if (.not.post_no_flow_mem_alloc) then

			call compute_mix_mean_vel

!		allocate (dissip(nx,my,mz),fluc(nx,my,mz))
!		dissip = 0d0
!		fluc   = 0d0

			if (near_particle_checked==.false.) then
				call near_particle_region
				call find_fluid_plus
				near_particle_checked = .true.
			endif

			sijsij_loc     = 0d0
			sijsij_far_loc = 0d0
			! generating thetaij
			if (I_AM_NODE_ZERO) write (*,*) "GENERATING SijSij"
			do i=1, nx
				if (I_AM_NODE_ZERO.and.mod(i,20)==0) write (*,*) " PLANE # = ", i
				do dim1=1, ndim
					do dim2=1, ndim
						call derivative_1st(i,dim1,dim2,ur1)
						call derivative_1st(i,dim2,dim1,ur2)

						do k=1, mz
							do j=1, my
								if (fluid_atijk(i,j,k)) then
									tmp1 = 2*vis*((ur1(j,k)+ur2(j,k))*half)**2
									sijsij_loc = sijsij_loc + tmp1
			
	!								dissip(i,j,k) = dissip(i,j,k) + tmp1
	!								fluc(i,j,k) = dot_product(ubcp(i,j,k,:)-ufmean(:),ubcp(i,j,k,:)-ufmean(:))/2

									if (far_field(i,j,k)) then
										sijsij_far_loc = sijsij_far_loc + tmp1
									endif
								else
	!								dissip(i,j,k) = -(vis*(umeanslip/dia_phys)**2)
	!								fluc(i,j,k)   = -(dot_product(ufmean(:),ufmean(:))/2)
								endif
							enddo
						enddo
					enddo
				enddo
			enddo
	!		dissip = dissip / (vis*(umeanslip/dia_phys)**2)
	!		fluc   = fluc   / (dot_product(ufmean(:)-usmean(:),ufmean(:)-usmean(:))/2)

#if PARALLEL
			sijsij = 0d0
			GLOBAL_DOUBLE_SUM(sijsij_loc,sijsij,1,decomp_group)

			sijsij_far = 0d0
			GLOBAL_DOUBLE_SUM(sijsij_far_loc,sijsij_far,1,decomp_group)
#else
			sijsij     = sijsij_loc
			sijsij_far = sijsij_far_loc
#endif
			sijsij_pure= sijsij     / count_fluid
			sijsij     = sijsij     / count_fluid / (vis*(mix_mean_slip_mag/dia_phys)**2)
			sijsij_far = sijsij_far / count_fluid / (vis*(mix_mean_slip_mag/dia_phys)**2)


			int_tke     = zero
			pres_part   = zero
			visc_part   = zero

			int_mom     = zero
			total_force = zero
			do m=1, nbody
				do idim=1, ndim
					tmp1 = ufmean(idim)-usmean(idim)
					pres_part   = pres_part + tmp1 * pres(m,idim)
					visc_part   = visc_part + tmp1 * visc(m,idim)
					int_tke     = int_tke   + tmp1 * (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)

					int_mom(idim)     = int_mom(idim)     + (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)
					total_force(idim) = total_force(idim) + (visc(m,idim)+pres(m,idim)-mpg(idim)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0)
				enddo
			enddo

			int_tke_fluc   = zero
			pres_part_fluc = zero
			visc_part_fluc = zero
			if (imove==1) then
				do m=1, nbody
					do idim=1, ndim
						tmp1 = velbdy(m,idim) - usmean(idim)

						pres_part_fluc = pres_part_fluc - pres(m,idim)                * tmp1
						visc_part_fluc = visc_part_fluc - visc(m,idim)                * tmp1
						int_tke_fluc   = int_tke_fluc   - (visc(m,idim)+pres(m,idim)) * tmp1
					enddo
				enddo
			endif

			!^^^ NORMALIZATION ^^^^^^
			int_mom(:)     = int_mom(:)     /phase_array(1)%npart/(3.d0*pi*vis*(mix_mean_slip_mag) * phase_array(iphs)%dia)
			total_force(:) = total_force(:) /phase_array(1)%npart/(3.d0*pi*vis*(mix_mean_slip_mag) * phase_array(iphs)%dia)

			int_mom_mag     = sqrt(dot_product(int_mom(:),int_mom(:)))
			total_force_mag = sqrt(dot_product(total_force(:),total_force(:)))

			pres_part = pres_part / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) / (vis*(mix_mean_slip_mag/dia_phys)**2)
			visc_part = visc_part / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) / (vis*(mix_mean_slip_mag/dia_phys)**2)

			int_tke_pure = int_tke   / (doml(1)*doml(2)*doml(3))
			int_tke      = int_tke   / (doml(1)*doml(2)*doml(3)) / (vis*(mix_mean_slip_mag/dia_phys)**2) !/ (1-maxvolfrac)

			pres_part_fluc = pres_part_fluc / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) / (vis*(mix_mean_slip_mag/dia_phys)**2)
			visc_part_fluc = visc_part_fluc / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) / (vis*(mix_mean_slip_mag/dia_phys)**2)
			int_tke_fluc   = int_tke_fluc   / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)) / (vis*(mix_mean_slip_mag/dia_phys)**2)

			if (I_AM_NODE_ZERO) then
				filename = trim(run_name)//"_sijsij.dat"
				open (unit=1, file=trim(filename), status="replace")

				if (.not.imove==1) then
					if (nphases==1) then
						write (1,"(4f8.4,6D15.7)") dbydx, lybyd, maxvolfrac, re,			&
											& sijsij, sijsij_far, sqrt(sijsij)/maxvolfrac,		&
											& int_tke, int_tke/nbody, int_tke_fluc
					else
						write (1,"(8f8.4,11D15.7)") dbydx, lybyd, maxvolfrac, (phase_array(iphs)%volfracg, iphs=1,nphases), 	&
											& (yalpha(iphs),iphs=1,nphases), re,										&
											& sijsij, sijsij_far, sijsij_far/sijsij,				&
											& (int_tke(iphs), iphs=1,nphases), mix_int_tke, 								&
											& (int_tke(iphs)/phase_array(iphs)%npart, iphs=1,nphases),						&
											& (int_tke_fluc(iphs), iphs=1,nphases), mix_int_tke_fluc
					endif
				else
					if (nphases==1) then
						write (1,"(6f8.4,6D15.7)") dbydx, lybyd, maxvolfrac, re, rhos/rhof, coeff_rest,	&
											& sijsij, sijsij_far, sijsij/sijsij_far,			&
											& int_tke, int_tke/nbody, int_tke_fluc
					else
						write (1,"(10f8.4,11D15.7)") dbydx, lybyd, maxvolfrac, (phase_array(iphs)%volfracg, iphs=1,nphases), &
											& (yalpha(iphs), iphs=1,nphases), re,										&
											& rhos, coeff_rest,														&
											& sijsij, sijsij_far, sijsij_far/sijsij,				&
											& (int_tke(iphs), iphs=1,nphases), mix_int_tke, 								&
											& (int_tke(iphs)/phase_array(iphs)%npart, iphs=1,nphases),						&
											& (int_tke_fluc(iphs), iphs=1,nphases), mix_int_tke_fluc
					endif
				endif
				close (1)
!
!				if (.not.imove==1) then
!					filename = trim(run_name)//"_force_comp.dat"
!					open (unit=1, file=trim(filename), status="replace")
!					if (nphases==1) then
!						write (1,"(4f8.4,2D15.7)") dbydx, lybyd, maxvolfrac, re,											&
!								& (int_mom_mag(iphs), iphs=1, nphases), (total_force_mag(iphs), iphs=1, nphases)
!					else
!						write (1,"(8f8.4,6D15.7)") dbydx, lybyd, maxvolfrac, (phase_array(iphs)%volfracg, iphs=1,nphases),	&
!								& (yalpha(iphs), iphs=1,nphases), re,											&
!								& (int_mom_mag(iphs), iphs=1, nphases), (total_force_mag(iphs), iphs=1, nphases), 		&
!								& mix_int_mom_mag, mix_total_force_mag
!					endif
!					close(1)
!				endif
			endif
		else
			if(I_AM_NODE_ZERO)then
				if (.not.imove==1) then
					if (nphases==1) then
						nvar  = 10
						nvar1 = 4
						nvar2 = 6
					else
						nvar  = 21
						nvar1 = 8
						nvar2 = 13
					endif
				else
					if (nphases==1) then
						nvar  = 12
						nvar1 = 6
						nvar2 = 6
					else
						nvar  = 23
						nvar1 = 10
						nvar2 = 13
					endif
				endif
				filename = "_sijsij.dat"
				call mis_average(nvar, nvar1, nvar2, filename, line)

				if (.not.imove==1) then
					if (nphases==1) then
						nvar  = 6
						nvar1 = 4
						nvar2 = 2
					else
						nvar  = 14
						nvar1 = 8
						nvar2 = 6
					endif
					filename = "_force_comp.dat"
					call mis_average(nvar, nvar1, nvar2, filename, line)
				endif
			endif 
		endif	
	end subroutine compute_sijsij
#endif

	subroutine mis_average(n, n1, n2, string, nline)
		implicit none

		integer, intent(in) :: n, n1, n2
		integer, intent(inout) :: nline
		character*50, intent(in) :: string

		integer :: imis,ivar, size, iline, io, nmis_tmp, nline_max=100000, nline_array(nmis)
		real(prcn), allocatable, dimension(:,:,:) :: var
		real(prcn), allocatable, dimension(:) :: avr_var, var_var
		real(prcn) :: confint
		character*50 filename, outform, outform_drag, zone
		character*2 tmp1, tmp2, tmp3, s1, s2, s3
		logical :: filexist
		real(prcn) :: fisol, ff, ffre, ft

		nline_array = 0

		call get_confin(nmis, confint)
		write (*,*) "CONFIDENCE = ", confint

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
		outform = trim(outform)//"f10.3,"
		outform = trim(outform)//trim(tmp2)
		outform = trim(outform)//"D15.7)"
		outform = "(2"//trim(outform)//")"

		outform = trim(outform)//"D15.7)"
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

		filename = "NMIS"//trim(string)
		open (unit=1,file=trim(filename),status="replace",action="write")

		allocate(avr_var(n),var_var(n))

		do iline=1, minval(nline_array(1:nmis_tmp))
			avr_var = 0d0
			var_var = 0d0

			do ivar=1, n
!				avr_var(ivar) = sum(var(:, iline, ivar))

				call calc_avr_var(nmis_tmp, var(1:nmis_tmp, iline, ivar), avr_var(ivar), var_var(ivar))
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

			write (1,trim(outform)) avr_var(:), var_var(:)
		enddo
		close (1)
10		if (allocated(var)) deallocate(var)
		if (allocated(var_var)) deallocate (avr_var, var_var)
	end subroutine mis_average


	subroutine calc_avr_var(nvar, var, avr, variance)
		implicit none
		integer, intent(in) :: nvar
		real(prcn), intent(in) :: var(nvar)
		real(prcn), intent(out) :: avr, variance

		integer :: ivar

		avr = sum(var(:))/nvar

		variance = zero
		do ivar=1, nvar
			variance = variance + (var(ivar)-avr)**2
		enddo
	end subroutine calc_avr_var




	subroutine derivative_1st(pos,dim1,dim2,duidxj)
		implicit none
		integer, intent(in)  :: pos, dim1, dim2
		real(prcn), intent(out) :: duidxj(my,mz)
	
		real(prcn) :: tmp1, tmp2
		integer :: i,j,k,im,ip
	
		i = pos
		ip=i+1
		im=i-1
	
#if !PARALLEL
		if(im<1)  im = nx
		if(ip>nx) ip = 1
#endif       
		if (dim2==1) then
			do k=1, mz
				do j=1, my
					duidxj(j,k)=(1.0/2.0/dx)* (ubcp(ip,j,k,dim1)-ubcp(im,j,k,dim1))
				enddo
			enddo
		else
			do k=1, mz
				do j=1, my2
					if (dim2==2) then
						uf3(j,k) = u(i,j,k,dim1) * wy(j)
					else
						uf3(j,k) = u(i,j,k,dim1) * wz(k)
					endif
				enddo
			enddo
			call ff2cr(uf3,ur3)
			duidxj = ur3
		endif
	end subroutine derivative_1st

	subroutine near_particle_region
		implicit none
	
		integer :: i, j, k, ii, jj, kk, m, idim, imin, imax, jmin, jmax, kmin, kmax
		real(prcn), dimension(ndim) :: xlr, xll
		integer, dimension(ndim) :: cor_min, cor_max
		real(prcn) :: dist

		if (allocated(far_field)) deallocate(far_field)

		allocate(far_field(0:nx+1,my,mz))
		do k = 1, mz
			do j = 1, my
				do i = 0, nx+1 
					far_field(i,j,k) = .true.
				end do
			end do
		end do

	    do m = 1, nbody
		  do idim = 1, ndim 
		     if(idim.eq.1) then 
		        xlr(idim) = xc(m,idim)  + radbdy(m) + foffset
		        xll(idim) = xc(m,idim)  - radbdy(m) + foffset
		     else 
		        xlr(idim) = xc(m,idim)  + radbdy(m)
		        xll(idim) = xc(m,idim)  - radbdy(m) 
		     end if
		  end do
	    
		  do idim = 1, ndim 
		     cor_min(idim)  = ceiling(xll(idim))
		     cor_max(idim) = floor(xlr(idim)) 
		  end do
	    
		  imin = cor_min(1)
		  imax = cor_max(1)
		  jmin = cor_min(2)
		  jmax = cor_max(2)
		  kmin = cor_min(3)
		  kmax = cor_max(3)
		  
		  do i = imin, imax 
		     ii = i
#if PARALLEL
		     if(.not.POINT_IN_PROC(ii))then

		        WEST_PERIODIC_POINT_IMAGE(i,ii)
		        EAST_PERIODIC_POINT_IMAGE(i,ii)
		        if(.not.POINT_IN_PROC(ii)) goto 555
		     end if
		     
#else
		     if(i.lt.1.and.intx_per) ii = mxf+i-1
		     if(i.gt.mxf-1.and.intx_per) ii = i-(mxf-1)
#endif
		     
		     LOCAL_INDEX(ii)

		     do j = jmin, jmax 
		        jj = j 
		        if(j.lt.1.and.inty_per) jj = my+j
		        if(j.gt.my.and.inty_per) jj = j-my
		        
		        do k = kmin, kmax 
		           
		           kk = k 
		           if(k.lt.1.and.intz_per) kk = mz+k
		           if(k.gt.mz.and.intz_per) kk = k-mz 
		           
		           dist = (i - (xc(m,1)+foffset))**two
		           dist = dist + (j - xc(m,2))**two + (k - xc(m,3))**two 
		           dist = dsqrt(dist)
		           IF((dist - radbdy(m)*near_particle_rad).LE.SMALL_NUMBER) THEN 
		              if(far_field(ii,jj,kk).eq..false.) then 
		                 PRINT*,'FLUID_ATIJK ALREADY FALSE AT I,J,K=&
		                      & ', ii,jj,kk, myid,m
		              ENDIF
		              far_field(ii,jj,kk)  = .false.
		           end IF
		        end do
		     end do
#if PARALLEL
555       continue
		     
#endif
		  end do
	    end do
	end subroutine near_particle_region


	subroutine find_fluid_plus
		implicit none
		integer :: i, j, k, ii, jj, kk, im, ip, jm, jp, km, kp
		logical :: bound

		if (allocated(fluid_atijkp)) deallocate(fluid_atijkp)

		allocate(fluid_atijkp(0:nx+1,my,mz))
		fluid_atijkp = fluid_atijk

		do k=1, mz
			do j=1, my
				do i=1, nx
					if (fluid_atijk(i,j,k)) then
						im = i-1
						ip = i+1
						if (i==1)  im = nx
						if (i==nx) ip = 1

						jm = j-1
						jp = j+1
						if (j==1)  jm = my
						if (j==my) jp = 1

						km = k-1
						kp = k+1
						if (k==1)  km = mz
						if (k==mz) kp = 1

						bound = .false.
						do kk=km, kp
							do jj=jm, jp
								do ii=im, ip
									if (.not.fluid_atijk(ii,jj,kk)) then
										bound = .true.
										exit
									endif
								enddo
								if (bound) exit
							enddo
							if (bound) exit
						enddo
						if (bound) fluid_atijkp(i,j,k) = .false.
					endif
				enddo
			enddo
		enddo
	end subroutine find_fluid_plus


	subroutine reynolds_stress_tensor
		implicit none
		real(prcn) :: spec_fluc(nphases), uiuj_s(nphases,ndim,ndim), tke_s(nphases), tke_mix_s
		real(prcn), dimension(ndim,ndim) :: uiuj_f, uiuj_far, bij_f, bij_s
		real(prcn) :: xi, eta1
		real(prcn) :: meanslip_energy, tke_far

		real(prcn) :: tmp_tensor(ndim,ndim), tmp_vec(ndim)
		integer :: dim1, dim2, dim3
		integer :: i, j, k, ip, im, idim, pstart, pend, m, iphs
		character*3 mis_str
		character*50 filename

		real(prcn), pointer, dimension(:,:,:,:) :: ur

		integer :: nvar, nvar1, nvar2
		logical :: filexist

		if (I_AM_NODE_ZERO) write (*,*) "IN REYNOLDS_STRESS_TENSOR"
		if (.not.post_no_flow_mem_alloc) then
			call compute_mix_mean_vel
			ur=>ubcp
			if (from_post) then
				if (near_particle_checked==.false.) then
					call near_particle_region
					call find_fluid_plus
					near_particle_checked = .true.
				endif
			endif

			!^^^^^^^ Computing RSM_F ^^^^^^^^^^^^^^
			uiuj_f = 0d0
			uiuj_far = 0d0
			do k=1, mz
				do j=1, my
					do i=1, nx
						if (fluid_atijk(i,j,k)) then
							do dim1=1, ndim
								do dim2=1, ndim
									uiuj_f(dim1,dim2) = uiuj_f(dim1,dim2) + (ubcp(i,j,k,dim1)-ufmean(dim1)) * (ubcp(i,j,k,dim2)-ufmean(dim2))
								enddo
							enddo

							if (from_post) then
								if (far_field(i,j,k)) then
									do dim1=1, ndim
										do dim2=1, ndim
											uiuj_far(dim1,dim2) = uiuj_far(dim1,dim2) + (ur(i,j,k,dim1)-ufmean(dim1)) * (ur(i,j,k,dim2)-ufmean(dim2))
										enddo
									enddo
								endif
							endif
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
							do dim2=1, ndim
								uiuj_s(iphs,dim1,dim2) = uiuj_s(iphs,dim1,dim2) + (velbdy(m,dim1)-phase_array(iphs)%mean_spec_vel(dim1)) * (velbdy(m,dim2)-phase_array(iphs)%mean_spec_vel(dim2))
							enddo
						end do
					end do
					uiuj_s(iphs,:,:) = uiuj_s(iphs,:,:)/real(phase_array(iphs)%npart,prcn)
					pstart = pend + 1
				end do

				do iphs = 1, nphases
					do idim=1, ndim
						tke_s(iphs) = tke_s(iphs) + uiuj_s(iphs,idim,idim)
					enddo
				end do

				tke_s(:) = tke_s(:)/2
			endif

			if (.not.allocated(tke_s_pure)) allocate(tke_s_pure(nphases))

			tke_s_pure = tke_s

			do iphs = 1, nphases
				tke_mix_s = tke_mix_s + tke_s(iphs)*phase_array(iphs)%volfracg/maxvolfrac
			end do
#if PARALLEL
			tmp_tensor = 0d0
			GLOBAL_DOUBLE_SUM(uiuj_f, tmp_tensor, 9, decomp_group)
			uiuj_f = tmp_tensor

			tmp_tensor = 0d0
			GLOBAL_DOUBLE_SUM(uiuj_far, tmp_tensor, 9, decomp_group)
			uiuj_far = tmp_tensor
#endif
 			meanslip_energy = half*mix_mean_slip_mag**2

			uiuj_f   = uiuj_f   / count_fluid
			uiuj_far = uiuj_far / count_fluid

			tke     = (uiuj_f(1,1)   + uiuj_f(2,2)   + uiuj_f(3,3))   / 2
			tke_far = (uiuj_far(1,1) + uiuj_far(2,2) + uiuj_far(3,3)) / 2

			tke_pure = tke

			uiuj_f = uiuj_f / meanslip_energy
			uiuj_far = uiuj_far / meanslip_energy

			tke     = tke     / meanslip_energy
			tke_far = tke_far / meanslip_energy

!			bij_f(:,:) = uiuj_f(:,:)/2.0/tke
!			do idim=1, ndim
!				bij_f(idim,idim) = bij_f(idim,idim)-1.0/3.0
!			enddo

			uiuj_s    = uiuj_s / meanslip_energy
			tke_s     = tke_s  / meanslip_energy
			tke_mix_s = tke_mix_s / meanslip_energy


			if(I_AM_NODE_ZERO)then
				call screen_separator(40,'K')

				do iphs = 1, nphases
					call screen_separator(20,'^')
					write (*,"(1A,1I3)") "SPECIES #", iphs
					write (*,"(1A,4D15.7)") "MEAN VELOCITY, MAG = ", phase_array(iphs)%mean_spec_vel(:), sqrt(dot_product(phase_array(iphs)%mean_spec_vel(:), phase_array(iphs)%mean_spec_vel(:)))
					write (*,"(1A)") "REYNOLDS STRESS"
					do idim=1, ndim
						write (*,"(3D15.7)") uiuj_s(iphs,idim,:)
					enddo
					call calc_anisotropy(uiuj_s, bij_s, xi, eta1) 
				
					write (*,*)
					write (*,"(1A,1D15.7)") "TKE = ", tke_s(iphs)
					call screen_separator(20,'-')
				enddo

				write (*,"(1A,1D15.7)") "MIX MEAN SLIP = ", mix_mean_slip_mag
				write (*,"(1A,1D15.7)") "MIX TKE SOLID = ", tke_mix_s

				call calc_anisotropy(uiuj_f, bij_f, xi, eta1) 
				write (*,'(A,3D15.7)') "TKE_F/MEAN_energy,XI, ETA = ", tke, xi, eta1
				write (*,'(A)') "ANISOTROPY TENSOR = "
				write (*,'(3D15.7)') ((bij_f(i,j), j=1,ndim), i=1,ndim)
				call screen_separator(40,'K')


				filename = trim(run_name)//"_tke"
				if (from_post) then
					filename=trim(filename)//"_post.dat"
					open (unit=1, file=trim(filename), status="replace", action="write")
					write (1,*) "zone"
				else
					filename=trim(filename)//".dat"
					if (irestart==0.and.first_pass) then
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
				endif


				if (from_post) then
					if (.not.imove==1) then
						if (nphases==1) then
							write (1,"(4F10.2,21D15.7)") dbydx, lybyd, maxvolfrac, re,	&
							& t/t_conv/maxvolfrac, tke, tke_far/tke, &
							& ((uiuj_f(i,j),i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim)
						else
							write (1,"(6f10.2,21D15.7)") maxvolfrac, &
							& (phase_array(iphs)%volfracg, iphs=1,nphases), &
							& (yalpha(iphs), iphs=1,nphases), re,	&
							& t/t_conv/maxvolfrac, tke, tke_far/tke, &
							& ((uiuj_f(i,j),i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim)
						endif
					else
						if (nphases==1) then
							write (1,"(6f10.2,31D15.7)") dbydx, lybyd, maxvolfrac, &
							& re, rhos, coeff_rest, &
							& t/t_conv/maxvolfrac, tke,tke_far/tke, &
							& ((uiuj_f(i,j),i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim), &
							& tke_s(1), &
							& ((uiuj_s(1,i,j),i=1,ndim),j=1,ndim)
						else
							write (1,"(8f10.2,42D15.7)") maxvolfrac,(phase_array(iphs)%volfracg, iphs=1,nphases),(yalpha(iphs), iphs=1,nphases), &
							& re, rhos, coeff_rest, &
							& t/t_conv/maxvolfrac, tke,tke_far/tke, &
							& ((uiuj_f(i,j),i=1,ndim),j=1,ndim), &
							& ((bij_f(i,j),i=1,ndim),j=1,ndim), &
							& tke_s(1), &
							& ((uiuj_s(1,i,j),i=1,ndim),j=1,ndim), &
							& tke_s(2), &
							& ((uiuj_s(2,i,j),i=1,ndim),j=1,ndim), &
							& tke_mix_s
						endif
					endif
				else
					write (1,"(11d15.7)") t/dia_phys*umeanslip/(1-maxvolfrac), tke, xi, eta1, bij_f(1,1), bij_f(2,2), bij_f(3,3), &
												& bij_f(1,2), bij_f(1,3), bij_f(2,3)
				endif
				close (1)
#if 0
				if (imove==1) then
					filename = trim(run_name)//"_tke_mix.dat"
					open (unit=1,file=trim(filename),status="replace")
					write (1,"(4f8.4,4D15.7)") maxvolfrac, re, rhos, coeff_rest, tke, tke_mix_s, (tke*(1-maxvolfrac)*rhof + tke_mix_s*maxvolfrac*rhos), (tke*(1-maxvolfrac)*rhof + tke_s*maxvolfrac*rhos) / ((1-maxvolfrac)*rhof + maxvolfrac*rhos)
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
				call mis_average(nvar, nvar1, nvar2, filename, line)
			endif 
		endif
	end subroutine reynolds_stress_tensor

	subroutine combine_history
		implicit none
		real(prcn), allocatable, dimension(:,:,:) :: var, var_fin
		real(prcn), allocatable, dimension(:) :: var_tmp
		real(prcn), allocatable, dimension(:,:) :: avr_var, var_var
		integer, allocatable, dimension(:) :: step

		real(prcn) :: time, confint
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

		call get_confin(nmis, confint)
		write (*,*) "CONFIDENCE = ", confint

		allocate(step(nmis))

		do ihist=1, 3
			if (ihist==1) then
				filename="tke_history.dat"
!				nvar = 24
!				des_var1 = 5
!				des_var2 = 5

				nvar = 5
				des_var1 = 1
				des_var2 = 2
			elseif (ihist==2) then
				filename="norm_drag.dat"
				if (nphases==1) then
					nvar = 11
				else
					nvar = 12 !4 + nphases + ndim * 2
				endif
				des_var1 = 1
				des_var2 = 2
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
!							var_tmp(des_var1) = var_tmp(des_var1) / (one-maxvolfrac)
						elseif (ihist==2) then
							if (nphases==1) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
							elseif (nphases==2) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
							endif
							var_tmp(des_var1) = var_tmp(des_var1) / (one-phiavg)
						elseif (ihist==3) then
							if (nphases==1) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
!								var_tmp(8) = var_tmp(8)**2 * 3 * rhos * maxvolfrac
							elseif (nphases==2) then
								read(1,*,iostat=io) var_tmp (1:nvar) ! var(imis,1:nvar,step(imis))
							endif
!							var_tmp(des_var1) = var_tmp(des_var1) / (one-maxvolfrac)
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

				write (*,*) "95% CONFIDENCE INTERVAL COMPUTED"
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


	subroutine correct_time
		implicit none
		integer :: nvar1, nvar2, imis, io1, io2, istep
		real(prcn), allocatable :: var1(:), var2(:)
		character*50 filename1, filename2, filename3
		logical :: filexist

		nvar1 = 23
		nvar2 = 11


		call screen_separator(30,'^')
		write (*,*) "IN CORRECT_TIME..."

		do imis=1, nmis
			if (imis==1) then
				filename1="MIS1_"//"tke_history.dat"
				filename2="MIS1_"//"norm_drag.dat"
				filename3="MIS1_"//"tke_history2.dat"
			elseif (imis==2) then
				filename1="MIS2_"//"tke_history.dat"
				filename2="MIS2_"//"norm_drag.dat"
				filename3="MIS2_"//"tke_history2.dat"
			elseif (imis==3) then
				filename1="MIS3_"//"tke_history.dat"
				filename2="MIS3_"//"norm_drag.dat"
				filename3="MIS3_"//"tke_history2.dat"
			elseif (imis==4) then
				filename1="MIS4_"//"tke_history.dat"
				filename2="MIS4_"//"norm_drag.dat"
				filename3="MIS4_"//"tke_history2.dat"
			elseif (imis==5) then
				filename1="MIS5_"//"tke_history.dat"
				filename2="MIS5_"//"norm_drag.dat"
				filename3="MIS5_"//"tke_history2.dat"
			elseif (imis==6) then
				filename1="MIS6_"//"tke_history.dat"
				filename2="MIS6_"//"norm_drag.dat"
				filename3="MIS6_"//"tke_history2.dat"
			elseif (imis==7) then
				filename1="MIS7_"//"tke_history.dat"
				filename2="MIS7_"//"norm_drag.dat"
				filename3="MIS7_"//"tke_history2.dat"
			elseif (imis==8) then
				filename1="MIS8_"//"tke_history.dat"
				filename2="MIS8_"//"norm_drag.dat"
				filename3="MIS8_"//"tke_history2.dat"
			endif
		
			inquire(file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				write (*,*) 'FILE "'//trim(filename1)//' DOES NOT EXIST, SKIPPING '//trim(filename1)
				goto 10
			endif

			write (*,*) "READING FROM "//trim(filename1)//" and "//trim(filename2)
			open (unit=1,file=trim(filename1),status="old",action="read")
			open (unit=2,file=trim(filename2),status="old",action="read")
			open (unit=3,file=trim(filename3),status="replace")

			if (.not.allocated(var1)) allocate(var1(nvar1), var2(nvar2))

			istep = 0
			do
				read (1,*,iostat=io1) var1(1:nvar1)
				read (2,*,iostat=io2) var2(1:nvar2)

				if (io1>0.or.io2>0) then
					write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN LINE ", istep+1
				elseif (io1==0.and.io2==0) then
					istep = istep+1
					write (3,"(4F8.4,20D15.7)") var1(1:4), var2(1), var1(5:23)
				elseif (io1<0.or.io2<0) then
					write (*,*) "END OF FILE "//trim(filename1)//" or "//trim(filename2)

					close(1)
					close(2)
					close(3)
					call screen_separator(30,'-')
					exit
				endif
			enddo
10			continue
		enddo
	end subroutine correct_time

#if 0
	SUBROUTINE GET_CONFIN(ndatain, confin)
		implicit none 
		INTEGER, INTENT(IN) :: ndatain 
		REAL(prcn), INTENT(OUT) ::  confin
		CHARACTER*50 :: table_name
		INTEGER :: tunit, ndata 
		LOGICAL :: filexist,isopen

		integer, parameter :: nmismax = 30, percmax = 11
		character*8 :: per_conf="95%"

		Type :: conf_interval
		 character*8 :: confper(percmax)
		 real(prcn) :: confi(percmax)
		END Type conf_interval

		TYPE(conf_interval) :: cis(nmismax)

		ndata = ndatain 
		!  WRITE(*,*)'IN GET CONFIN, NDATA = ', ndata
		table_name = "ci_lookup.dat"
		INQUIRE(FILE=trim(table_name),EXIST=filexist,OPENED=isopen)
		if(filexist)then
			tunit = getnewunit(minunitno,maxunitno)
			OPEN(unit=tunit,FILE=table_name,form='formatted')
			CALL initialize_lookup_table(tunit)
			if(ndata.gt.nmismax)ndata=nmismax
			if(ndata.gt.1)then
			confin = lookup_entry(ndata-1,per_conf)

			!WRITE(*,'(2(A,2x,g17.8))') 'DATA PTS = ', real(NDATA),' CONFIN = ', CONFIN
			else
			!WRITE(*,'(A,2x,g17.8)') 'ONLY ONE MIS: CONFIN = ', CONFIN
			confin = zero
			end if
			close(tunit,status='keep')

		else

			write (*,*) trim(table_name)//" DOES NOT EXIST. CHECK THE FILE...."
			stop
			confin = 1.812
		end if

	contains
		SUBROUTINE initialize_lookup_table(tunit)
			IMPLICIT NONE	


			CHARACTER*8 :: ciperc(percmax)
			INTEGER, Intent(in) :: tunit
			INTEGER :: iperc,isamples

			READ(tunit,*)(ciperc(iperc),iperc=1,percmax)

			do isamples=1,nmismax
				cis(isamples)%confper(1:percmax) = ciperc(1:percmax)
			end do
			do isamples=1,nmismax
				READ(tunit,*)(cis(isamples)%confi(iperc),iperc=1,percmax)
			end do

		END SUBROUTINE initialize_lookup_table

		RECURSIVE function lookup_entry(nsamples,entry) result(confis)
			IMPLICIT NONE
			INTEGER, Intent(in) :: nsamples
			CHARACTER*8, Intent(inout) :: entry
			REAL(prcn) :: confis
			Integer :: iperc

			do iperc = 1,percmax
				if(entry.eq.(cis(nsamples)%confper(iperc)))then
					confis = cis(nsamples)%confi(iperc)
					exit
				end if
			end do

			if(iperc.gt.percmax)then
			!       PRINT*,'ENTRY ', entry,' NOT FOUND. SETTING 95% CONFIDENCE INTERVAL'
				entry = TRIM('95%')
				confis = lookup_entry(nsamples,entry)
				RETURN
			end if
!    PRINT*,'CONFIDENCE INTERVAL FOR ', TRIM(entry),' CONFIDENCE AND ',nsamples, ' DEGREES OF FREEDOM IS ', confis 
			RETURN
		END function lookup_entry
	end SUBROUTINE GET_CONFIN
#endif
	subroutine cluster_characteristics
		implicit none

		character*50 filename1, filename2, str_tmp
		integer :: idim, jdim, ibody
		
		integer :: i, j, k, ibin, nvar, nvar1, nvar2
		real(prcn) :: dr, rmax, max_overlap, tmp1, tmp2, norm_force
		real(prcn), allocatable :: gofr(:), force_aggr(:,:)

		integer :: icluster
		logical, allocatable :: contact(:,:)


		if (.not.I_AM_NODE_ZERO) return

		if (.not. post_no_flow_mem_alloc) then
			call screen_separator(30,'^')
			write (*,*) "IN CLUSTER_CHARACTERISTICS..."

			call compute_mix_mean_vel
			allocate(contact(nbody,nbody), force_aggr(nbody,ndim))
			contact = .false.

			if (allocated(rad_bin)) deallocate(rad_bin)
			allocate(rad_bin(1:nbins))
			allocate(gofr(1:nbins))
			rad_bin = zero
			gofr = zero

			rmax  = lybyd / 2
			dr = rmax / nbins
			do ibin=1, nbins
				rad_bin(ibin) = (ibin-0.5) * dr
			enddo

			call calculate_gofr_homog(nbody, xc(:,:)/dbydx, contact, my, mx+1, nbins, .true., gofr(1:nbins), rad_bin(1:nbins), max_overlap)

			force_aggr(:,:) = force(:,:)
	!		if (small_max_overlap>max_overlap) write (*,"(1a,1d15.7)") "max overlap = ", small_max_overlap

			filename2 = trim(run_name)//"_cluster"
			call find_clusters(xc/dbydx, contact, nbody, lybyd, filename2, .true.)


			filename1 = trim(run_name)//"_config.dat"
			open (unit=1, file=trim(filename1), status="replace", action="write")
			write (1,*) "zone"

			norm_force = 3.d0*pi*vis*mix_mean_slip_mag*(1-maxvolfrac)*dia_phys
			do ibody=1, nbody
				tmp1 = sqrt(dot_product(force(ibody,:), force(ibody,:)))
				tmp2 = sqrt(dot_product(force_aggr(ibody,:), force_aggr(ibody,:)))


				write (1,"(4d15.7,1i,2d15.7)") xc(ibody,1:ndim)/dbydx, radbdy(ibody)/dbydx, aggregate_num(ibody), tmp1/norm_force, tmp2/norm_force
			enddo
			close (1)

			call compute_force

			call compute_force_distribution_per_particle


			write (*,*) "LEAVING CLUSTER_CHARACTERISTICS..."
			call screen_separator(30,'-')
		else
			nvar  = 6
			nvar1 = 1
			nvar2 = 5
!			line = nbins

			filename1 = "_force_dist.dat"
			call mis_average(nvar, nvar1, nvar2, filename1, line)


			nvar  = 3
			nvar1 = 2
			nvar2 = 1
!			line = nbins

			filename1 = "_cluster_dist.dat"
			call mis_average(nvar, nvar1, nvar2, filename1, line)

			nvar  = 4
			nvar1 = 3
			nvar2 = 1
!			line = nbins

			filename1 = "_force_pdf.dat"
			call mis_average(nvar, nvar1, nvar2, filename1, line)


			nvar  = 4
			nvar1 = 3
			nvar2 = 1
!			line = nbins

			filename1 = "_force_pdf_conditinal.dat"
			call mis_average(nvar, nvar1, nvar2, filename1, line)
		endif

	contains

		subroutine compute_force
			implicit none
			integer :: i, j, idim, max_num, ibin, max_part
			real(prcn), allocatable :: force_bin(:), count_bin(:)
			real(prcn), allocatable :: force_bin_val(:,:), force_bin_mag(:), equi_part_dia1(:), equi_part_force1(:), equi_part_dia2(:), equi_part_force2(:)
			real(prcn) :: total_force(ndim), total_force_mag, norm_force, tmp1
			character*50 filename1


			max_part = maxval(cluster(:)%num)
!			max_num = max_part
			max_num = max_part/5+1

			allocate(force_bin(max_num), count_bin(max_num))
			allocate(force_bin_val(max_num,ndim), force_bin_mag(max_num))
			allocate(equi_part_dia1(max_num), equi_part_force1(max_num))
			allocate(equi_part_dia2(max_num), equi_part_force2(max_num))

			force_bin = 0
			force_bin_val = zero
			force_bin_mag = zero
			equi_part_dia1 = zero
			equi_part_dia2 = zero
			equi_part_force1 = zero
			equi_part_force2 = zero
			total_force = zero
			count_bin = zero

			do i=1, ncluster
!				ibin = cluster(i)%num
				ibin = cluster(i)%num/5+1

				if (ibin>max_num) ibin = max_num

				count_bin(ibin) = count_bin(ibin)+1
				force_bin(ibin) = force_bin(ibin) + cluster(i)%num

				do j=1, cluster(i)%num
					force_bin_val(ibin,:) = force_bin_val(ibin,:) + force(cluster(i)%list(j),:)

					total_force(:) = total_force(:) + force(cluster(i)%list(j),:)
				enddo
			enddo

			do ibin=1, max_num
				force_bin_mag(ibin) = sqrt(dot_product(force_bin_val(ibin,:), force_bin_val(ibin,:)))
			enddo
			total_force_mag = sqrt(dot_product(total_force(:), total_force(:)))


			do ibin=1, max_num
				equi_part_dia1(ibin) = dia_phys * (force_bin(ibin))**(1./3.)
				equi_part_dia2(ibin) = dia_phys * (force_bin(ibin)/count_bin(ibin))**(1./3.)
				call compute_ibm_drag(maxvolfrac, re*equi_part_dia1(ibin)/dia_phys, equi_part_force1(ibin))
				call compute_ibm_drag(maxvolfrac, re*equi_part_dia2(ibin)/dia_phys, equi_part_force2(ibin))
			enddo

			norm_force = 3.d0*pi*dia_phys*vis*(1-maxvolfrac)*mix_mean_slip_mag
			filename1 = trim(run_name)//"_force_dist.dat"
			open (unit=1, file=trim(filename1), status="replace", action="write")
			write (1,*) "zone"

			write (*,"(1a,1d15.7)") "FORCE CHECK = ", sum(force_bin_mag(:))/total_force_mag
			do ibin=1, max_num
				if (force_bin(ibin)>0) then
!					write (1,"(1f10.2, 4d15.7)") ibin, force_bin_mag(ibin)/force_bin(ibin)/norm_force, force_bin_mag(ibin)/norm_force, force_bin_mag(ibin)/total_force_mag, equi_part_force(ibin)
					write (1,"(1f10.2, 5d15.7)") 5*(ibin-.5), force_bin_mag(ibin)/force_bin(ibin)/norm_force, force_bin_mag(ibin)/norm_force, force_bin_mag(ibin)/total_force_mag, equi_part_force1(ibin), equi_part_force2(ibin)
				else
					write (1,"(1f10.2, 5d15.7)") 5*(ibin-.5), zero, zero, zero, zero, zero
				endif
			enddo
			close(1)



			filename1 = trim(run_name)//"_force_scatter.dat"
			open (unit=1, file=trim(filename1), status="replace", action="write")
			write (1,*) "zone"

			do i=1, ncluster
!				tmp1 = sqrt(dot_product())
!				write (1,"()") i, cluster(i)%num, 



			enddo








		end subroutine compute_force

		subroutine compute_force_distribution_per_particle
			implicit none

			real(prcn) :: norm_force, max_force, delta_force, max_part, delta_part
			real(prcn) :: pres_part(nbody,ndim), pres_part_mag(nbody), visc_part_mag(nbody), total_mag(nbody)
			integer :: ibody, ibin, jbin, bins, cluster_bins
			real(prcn), allocatable, dimension(:) :: pres_bin, visc_bin, total_bin
			real(prcn), allocatable, dimension(:,:) :: pres_bin_cluster, visc_bin_cluster, total_bin_cluster
			character*80 filename1

			pres_part = zero
			do ibody=1, nbody
				pres_part(ibody,:) = pres(ibody,:) - mpg(:)*pi*((two*radbdy(ibody)*dx)**3.d0)/6.d0
			enddo
    
			norm_force = 3.d0*pi*dia_phys*vis*(1-maxvolfrac)*mix_mean_slip_mag

			do ibody=1, nbody
				pres_part_mag(ibody) = sqrt(dot_product(pres_part(ibody,:),pres_part(ibody,:))) / norm_force
				visc_part_mag(ibody) = sqrt(dot_product(visc(ibody,:),visc(ibody,:)))           / norm_force
				total_mag(ibody)     = sqrt(dot_product(force(ibody,:),force(ibody,:)))         / norm_force
			enddo

			filename1 = trim(run_name)//"_force_per_part_scatter.dat"
			open (unit=1, file=trim(filename1), status="replace")
			write (1,*) "zone"

			do ibody = 1, nbody
				write (1,"(1i,4d15.7)") ibody, aggregate_num(ibody), pres_part_mag(ibody), visc_part_mag(ibody), total_mag(ibody)
			enddo
			close(1)

			bins = 10

			max_force = maxval(total_mag(:))
			delta_force = max_force/bins

			allocate(pres_bin(bins), visc_bin(bins), total_bin(bins))
			pres_bin = zero
			visc_bin = zero
			total_bin = zero

			do ibody=1, nbody
				ibin = int(pres_part_mag(ibody)/delta_force)+1
				if (ibin>bins) ibin = bins
				pres_bin(ibin) = pres_bin(ibin)+1

				ibin = int(visc_part_mag(ibody)/delta_force)+1
				if (ibin>bins) ibin = bins
				visc_bin(ibin) = visc_bin(ibin)+1

				ibin = int(total_mag(ibody)/delta_force)+1
				if (ibin>bins) ibin = bins
				total_bin(ibin) = total_bin(ibin)+1
			enddo

			if (sum(pres_bin(:))>small_number) pres_bin(:) = pres_bin(:)/sum(pres_bin(:))
			if (sum(visc_bin(:))>small_number) visc_bin(:) = visc_bin(:)/sum(visc_bin(:))
			if (sum(total_bin(:))>small_number) total_bin(:) = total_bin(:)/sum(total_bin(:))

			filename1 = trim(run_name)//"_force_pdf.dat"
			open (unit=1, file=trim(filename1), status="replace")
			write (1,*) "zone"
			do ibin=1, bins
				write (1,"(4d15.7)") (ibin-.5)/bins * max_force, pres_bin(ibin)/delta_force, visc_bin(ibin)/delta_force, total_bin(ibin)/delta_force
			enddo
			close(1)


			cluster_bins = 6
!			max_part = maxval(cluster(:)%num)
			max_part = 30
			delta_part = max_part/cluster_bins

			allocate(pres_bin_cluster(cluster_bins,bins), visc_bin_cluster(cluster_bins,bins), total_bin_cluster(cluster_bins,bins))
			pres_bin_cluster = zero
			visc_bin_cluster = zero
			total_bin_cluster = zero


			if (ncluster>0) then
				do ibody=1, nbody
					ibin = int(aggregate_num(ibody)/delta_part)+1
					if (ibin>cluster_bins) ibin = cluster_bins

					jbin = int(pres_part_mag(ibody)/delta_force)+1
					if (jbin>bins) jbin = bins
					pres_bin_cluster(ibin,jbin) = pres_bin_cluster(ibin,jbin)+1

					jbin = int(visc_part_mag(ibody)/delta_force)+1
					if (jbin>bins) jbin = bins
					visc_bin_cluster(ibin,jbin) = visc_bin_cluster(ibin,jbin)+1

					jbin = int(total_mag(ibody)/delta_force)+1
					if (jbin>bins) jbin = bins
					total_bin_cluster(ibin,jbin) = total_bin_cluster(ibin,jbin)+1
				enddo

				do ibin=1, cluster_bins
					if (sum(pres_bin_cluster(ibin,:))>small_number) &
						pres_bin_cluster(ibin,:) = pres_bin_cluster(ibin,:) / sum(pres_bin_cluster(ibin,:))
					if (sum(visc_bin_cluster(ibin,:))>small_number) &
						visc_bin_cluster(ibin,:) = visc_bin_cluster(ibin,:) / sum(visc_bin_cluster(ibin,:))
					if (sum(total_bin_cluster(ibin,:))>small_number) &
						total_bin_cluster(ibin,:) = total_bin_cluster(ibin,:) / sum(total_bin_cluster(ibin,:))
				enddo
				
				filename1 = trim(run_name)//"_force_pdf_conditinal.dat"
				open (unit=1, file=trim(filename1), status="replace")
				write (1,*) "zone"
				do ibin=1, cluster_bins
					do jbin=1, bins
						write (1,"(4d15.7)") (jbin-.5)/bins * max_force, pres_bin_cluster(ibin,jbin)/delta_force, visc_bin_cluster(ibin,jbin)/delta_force, total_bin_cluster(ibin,jbin)/delta_force
					enddo
				enddo
				close(1)
			endif

		end subroutine compute_force_distribution_per_particle

	end subroutine cluster_characteristics

	subroutine output_particle_config1(xc, color1, rad, nbody, filename1)
		implicit none
		integer, intent(in) :: nbody
		real(prcn), intent(in) :: xc(nbody, ndim), rad(nbody)
		integer, intent(in) :: color1(nbody)
!		logical, intent(in) :: first_pass
		character*50, intent(in) ::  filename1
		integer :: i, j, k, m

		if (first_pass) then
			open (unit=1, file=trim(filename1), status="replace")
		else
			open (unit=1, file=trim(filename1), status="old", position="append")
		endif

		write (1,*) "zone"
		write (1,"(4d15.7,1i)") ((xc(i,:), rad(i), color1(i)), i=1, nbody)
		close (1)
	end subroutine output_particle_config1

	subroutine output_particle_config2(xc, color1, color2, rad, nbody, filename1)
		implicit none
		integer, intent(in) :: nbody
		real(prcn), intent(in) :: xc(nbody, ndim), rad(nbody)
		integer, intent(in) :: color1(nbody), color2(nbody)
!		logical, intent(in) :: first_pass
		character*50, intent(in) ::  filename1
		integer :: i, j, k, m

		if (first_pass) then
			open (unit=1, file=trim(filename1), status="replace")
		else
			open (unit=1, file=trim(filename1), status="old", position="append")
		endif

		write (1,*) "zone"
		write (1,"(4d15.7,2i)") ((xc(i,:), rad(i), color1(i), color2(i)), i=1, nbody)
		close (1)
	end subroutine output_particle_config2
end module mypost_process

