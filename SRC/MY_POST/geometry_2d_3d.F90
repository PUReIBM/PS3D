module geometry_2d_3d
#include "../FLO/ibm.h"
	use global_data
	use nlmainarrays , only : ubcp, pbcp
	use fftw3_interface
	use constants
	use nlarrays
	use randomno
	use general_funcs
	use string_funcs
!	use boundary_condition
	use dependent_functions
	use bcsetarrays
!	use epsarray
!	use string_funcs
!	use init_turb
	use mypost_process
	use postproc_funcs
	use clustermod
	implicit none

	real(prcn), allocatable :: xc_conf(:,:), xc_conf_old(:,:), rad(:)

!	real(prcn), allocatable :: rad_bin(:)
	real(prcn), allocatable :: gofr_2d(:), gofr_var(:), gofr_2d_main(:), gofr_3d(:), gofr_3d_main(:)
	real(prcn), allocatable :: cofr_2d(:), cofr_var(:), cofr_2d_main(:), cofr_3d(:), cofr_3d_main(:)
	real(prcn), allocatable :: cofrf_2d(:), cofrf_var(:)

	real(prcn), allocatable :: cofr(:,:), cofrf(:,:), gofr(:,:)

	real(prcn), allocatable :: gofr_2d_coh(:,:)

	logical, allocatable :: center_changed(:)

	real(prcn) :: temperature

!	logical :: first_pass = .true.

contains
	subroutine rand_out(rn)
		implicit none

		real(prcn), intent(out) :: rn

		integer, save :: ran_max=10000

		real(prcn), save, allocatable :: ran_num(:)
		integer, save :: ran_count

		if (ran_count==ran_max) deallocate(ran_num)

		if (.not.allocated(ran_num)) then
			ran_count = 0
			allocate(ran_num(ran_max))
			call uni_dist(ran_num)
		endif

		ran_count = ran_count+1
		rn = ran_num(ran_count)
	end subroutine rand_out

	subroutine cohesive_particles
		implicit none

		character*50 filename1, filename2, filename3, filename4, str_tmp, filename5, filename6
		logical :: filexist
		real(prcn), allocatable :: part_tmp(:,:,:), part(:,:,:)
		real(prcn), allocatable :: gofr_coh(:,:), gofr_coh_avg(:), gofr_coh_var(:), e_par(:), e_per(:)
		real(prcn) :: xmin(ndim), xmax(ndim), vec1(ndim), vec2(ndim), vec3(ndim), rotation_tensor(ndim,ndim), r, confint
		real(prcn) :: point1(ndim), point2(ndim), point3(ndim)
		integer :: total_zone_max=50000000, total_zone, nzone, io, idim, jdim, ibody, jbody, kbody, izone, ibin, dim(2), dim1, dim2
		
		integer :: nsmall_box, iconfig, icluster, i, j, k, pcount, nhist
		real(prcn) :: volfrac_avg, volfrac_var, rmax, max_overlap, small_max_overlap, surface(ndim), cord(ndim), dist(2), sep, alpha
		real(prcn), allocatable :: config_tmp(:,:), vel_tmp(:,:), rad_tmp(:), gofr_small(:,:), gofr_small_avg(:), gofr_small_var(:), hist(:), bin_hist(:), cluster2d(:,:)

		logical, allocatable :: contact(:,:)

		integer :: turn, s_size
		character*4 turn_string


		type :: particle_type
			real(prcn) :: xc(ndim), rad
		end type particle_type

		type :: config_type
			real(prcn) :: orientation(ndim,ndim), center(ndim), volfrac, node(8,ndim)
			integer :: nbody

			type(particle_type), allocatable :: part(:)
		end type config_type
		type(config_type), allocatable :: small_config(:)

		logical :: inbox


		call screen_separator(30,'^')
		write (*,*) "IN COHISIVE_PARTS..."


		filename1 = "OUTPUT/output-001.dump"
		!filename1 = "prof19.dat"
!		filename2 = "cohesive_part_confing.dat"
		filename3 = "cohesive_gofr_energy.dat"
!		filename4 = "cohesive_cluster"
		filename5 = "cohesive_grant.dat"

		open (unit=3, file=trim(filename3), status="replace", action="write")
		open (unit=5, file=trim(filename5), status="replace", action="write")

		inquire(file=trim(filename1), exist=filexist)
		if (.not.filexist) then
			write (*,*) trim(filename1)//' DOES NOT EXIST. STOP'
			stop
		endif

		write (*,*) trim(filename1)//" IS AVAILABLE, ALLOCATING VARIABLES..."
		open (unit=1,file=trim(filename1),status="old",action="read")


		if (allocated(rad_bin)) deallocate(rad_bin)
		allocate(rad_bin(nbins))
		rad_bin = zero

		nzone = 0
		total_zone = 0 
		first_pass = .true.
		max_overlap = zero
		do

			if (mod(nzone,100)==0) then
				turn = nzone/100

				call to_string(turn, turn_string, s_size)

				do idim=1, 4-s_size
					turn_string="0"//trim(turn_string)
				enddo
				filename2 = "cohesive_part_confing"//trim(turn_string)//".dat"
				filename4 = "cohesive_cluster"//trim(turn_string)
				first_pass = .true.
			endif

			read (1,*,iostat=io) str_tmp
			if (io<0) then
				write (*,*) "END OF FILE ", trim(filename1)
				write (*,*) "NUMBER OF ZONES = ", nzone
				write (*,*) "NUMBER of PARTICLES = ", nbody
				exit
			endif
			total_zone = total_zone + 1
!			if (mod(total_zone,10)==1)
			write (*,*) "TOTAL_ZONE =", total_zone

			read (1,*) t
			read (1,*) str_tmp
			read (1,*) nbody
			read (1,*) str_tmp

			do idim=1, ndim
				read (1,*) xmin(idim), xmax(idim)
				lybyd = xmax(idim)-xmin(idim)
				rmax = lybyd/2
			enddo
			nx = lybyd !*dbydx
			my = nx
			mz = nx
			dia_phys = one

			read (1,*) str_tmp

			if (.not.allocated(config_tmp)) then
!				allocate(part_tmp(nbody, ndim+1))
				allocate(gofr_coh(1,nbins))
				allocate(e_par(nbins), e_per(nbins))

				allocate(contact(nbody,nbody))
				allocate(config_tmp(nbody,ndim), vel_tmp(nbody,ndim), rad_tmp(nbody))
			endif

			do ibody=1, nbody
				read (1,*) config_tmp(ibody,1:ndim), vel_tmp(ibody,1:ndim), rad_tmp(ibody)
			enddo

			if (mod(total_zone,skip_num)==0) then
				nzone = nzone + 1
				write (*,*) "NZONE = ", nzone

				contact = .false.
				call calculate_gofr_energy(nbody, config_tmp, vel_tmp, contact, my, mx+1, nbins, .true., gofr_coh(1,1:nbins), e_par, e_per, rad_bin, max_overlap)

				if (.not.allocated(grant)) allocate(grant(1))
				grant = zero
				do ibody = 1, nbody
					grant(1) = grant(1) + dot_product(vel_tmp(ibody,:),vel_tmp(ibody,:))
				end do
				grant(1) = grant(1)/(three*nbody)

				write (3,*) "zone"
				write (3,"(4d15.7)") ((rad_bin(ibin), gofr_coh(1,ibin), e_par(ibin)/(2*grant), e_per(ibin)/(2*grant)), ibin=1, nbins)

				write (5,"(2e15.7)") t, grant(1)

				call find_clusters(config_tmp, contact, nbody, lybyd, filename4, .true.)
				call output_particle_config1(config_tmp, aggregate_num(1:nbody), rad_tmp, nbody, filename2)
				write (*,"(1a,1d15.7)") "MAX OVERLAP = ", max_overlap
				first_pass = .false.
			endif
			!if (nzone>=10) exit
		enddo
		close(1)
		close(3)
		close(5)





#if 0

		filename1 = "OUTPUT/PosVelData-001.dump"
		!filename1 = "prof19.dat"

		inquire(file=trim(filename1), exist=filexist)
		if (.not.filexist) then
			write (*,*) trim(filename1)//' DOES NOT EXIST. STOP'
			stop
		endif

		write (*,*) trim(filename1)//" IS AVAILABLE, ALLOCATING VARIABLES..."
		open (unit=1,file=trim(filename1),status="old",action="read")

		nzone = 0
		total_zone = 0 
		do
			read (1,*,iostat=io) str_tmp
			if (io<0) then
				write (*,*) "END OF FILE ", trim(filename1)
				write (*,*) "NUMBER OF ZONES = ", nzone
				write (*,*) "NUMBER of PARTICLES = ", nbody
				exit
			endif
			total_zone = total_zone + 1
			if (mod(izone,10)==1) write (*,*) "ZONE ", nzone

			read (1,*) str_tmp
			read (1,*) str_tmp
			read (1,*) nbody
			read (1,*) str_tmp

			do idim=1, ndim
				read (1,*) xmin(idim), xmax(idim)
			enddo
			read (1,*) str_tmp

			if (.not.allocated(part_tmp)) allocate(part_tmp(total_zone_max/skip_num+1, nbody, ndim+1))

			do ibody=1, nbody
				read (1,*) vec1(:), vec2(:), r

				if (mod(izone,skip_num)==1) then
					nzone = nzone + 1

					part_tmp(nzone,ibody,1:ndim) = vec1(1:ndim) !+ 1
					part_tmp(nzone,ibody,4) = r
				endif
			enddo
		enddo
		close(1)

		allocate(part(nzone, nbody, ndim+1))
		do izone=1, nzone
			do ibody=1, nbody
				part(izone, ibody,:) = part_tmp(izone, ibody,:)
			enddo
		enddo
		deallocate(part_tmp)

		xmin = xmin !+ 1
		xmax = xmax !+ 1

		lybyd = xmax(1)-xmin(1)

		nx = lybyd !*dbydx
		my = nx
		mz = nx
		dia_phys = one

		rmax = lybyd/2

		write (*,*) "LYBYD = ", lybyd
		write (*,*) "NX = ", nx

!		do izone=1, nzone
!			do ibody=1, nbody
!				part(izone,ibody,:) = part(izone,ibody,:) * dbydx
!			enddo
!		enddo

!goto 20

		allocate(gofr_coh(nzone,nbins), gofr_coh_avg(nbins), gofr_coh_var(nbins))

		allocate(contact(nbody,nbody))
		allocate(config_tmp(nbody,ndim))
		config_tmp = zero

		gofr_coh = zero
		gofr_coh_avg = zero
		gofr_coh_var = zero

		if (allocated(rad_bin)) deallocate(rad_bin)

		allocate(rad_bin(nbins))
		rad_bin = zero

		filename2 = "OUTPUT/cohesive_part_confing.dat"
		filename3 = "OUTPUT/cohesive_gofr.dat"
		filename4 = "OUTPUT/cohesive_cluster"

		open (unit=3, file=trim(filename3), status="replace", action="write")

		write (*,*) "CALCULATING g(r) ...."

		first_pass = .true.
		max_overlap = zero
		do izone=1, nzone
!			if (mod(izone,skip_num)==1) then
				write (*,*) "ZONE = ", izone, "/", nzone

				contact = .false.

				config_tmp = part(izone,1:nbody,1:ndim)
				call calculate_gofr_homog(nbody, config_tmp, contact, my, mx+1, nbins, .true., gofr_coh(izone,1:nbins), rad_bin, max_overlap)
				call find_clusters(config_tmp, contact, nbody, lybyd, filename4, .true.)

				call output_particle_config1(part(izone,1:nbody,1:ndim), aggregate_num(1:nbody), part(izone,1:nbody,ndim+1), nbody, filename2)
				first_pass = .false.

				write (3,*) "zone"
				write (3,"(2d15.7)") ((rad_bin(ibin), gofr_coh(izone,ibin)), ibin=1,nbins)
!			endif
		enddo
		close(3)

		write (*,"(1a,1d15.7)") "MAX OVERLAP = ", max_overlap

		deallocate(gofr_coh, gofr_coh_avg, gofr_coh_var, contact, aggregate_num, config_tmp, rad_bin)
#endif


return


		if (nsmall_length> 1.and. (cluster_part_min==0.and.cluster_part_max==0)  ) then
			nsmall_box = nsmall_length**3
			lybyd_small = lybyd / nsmall_length

			allocate(small_config(nsmall_box))
			iconfig = 0
			do k=1, nsmall_length
				do j=1, nsmall_length
					do i=1, nsmall_length
						iconfig = iconfig+1

						small_config(iconfig)%center(1) = (i-0.5) * lybyd_small
						small_config(iconfig)%center(2) = (j-0.5) * lybyd_small
						small_config(iconfig)%center(3) = (k-0.5) * lybyd_small

						small_config(iconfig)%orientation = zero
						do idim=1, ndim
							small_config(iconfig)%orientation(idim,idim) = one
						enddo
					enddo
				enddo
			enddo
		elseif (cluster_part_min>0.and.cluster_part_max>0) then
			if (lybyd_small==one) then
				write (*,*) "LYBYD_SMALL NOT DEFINED"
				stop
			endif

			nsmall_box = 0
			do icluster=1, ncluster
				if (cluster_part_min<=cluster(icluster)%num .and. cluster(icluster)%num<=cluster_part_max) nsmall_box = nsmall_box+1
			enddo

			allocate(small_config(nsmall_box))

			iconfig = 0
			do icluster=1, ncluster
				if (cluster_part_min<=cluster(icluster)%num .and. cluster(icluster)%num<=cluster_part_max) then
					iconfig = iconfig+1

!					allocate(small_config(iconfig)%node(8,3))
!					small_config(nsmall_box)%node(8,3) = zero

					if (.not.rotated) then
!						allocate(small_config(nsmall_box)%start(ndim), small_config(nsmall_box)%end(ndim))

						small_config(iconfig)%orientation = zero
						do idim=1, ndim
							small_config(iconfig)%orientation(idim,idim) = one
						enddo

!						small_config(nsmall_box)%start(:) = cluster(icluster)%center(:)-lybyd_small/2
!						small_config(nsmall_box)%end(:)   = cluster(icluster)%center(:)+lybyd_small/2

					else
						do idim=1, ndim
							small_config(iconfig)%orientation(:,idim) = cluster(icluster)%eigen_vec(:,ndim-(idim-1))


						enddo
					endif

					small_config(iconfig)%center(:) = cluster(icluster)%center(:)

					point1(:) = cluster(icluster)%center(:) - small_config(iconfig)%orientation(:,1) * lybyd_small/2
					point2(:) = point1(:) - small_config(iconfig)%orientation(:,2) * lybyd_small/2
					point3(:) = point2(:) - small_config(iconfig)%orientation(:,3) * lybyd_small/2

					small_config(iconfig)%node(1,:) = point3(:)
					small_config(iconfig)%node(2,:) = small_config(iconfig)%node(1,:) + &
																	small_config(iconfig)%orientation(:,1) * lybyd_small
					small_config(iconfig)%node(3,:) = small_config(iconfig)%node(2,:) + &
																	small_config(iconfig)%orientation(:,2) * lybyd_small
					small_config(iconfig)%node(4,:) = small_config(iconfig)%node(1,:) + &
																	small_config(iconfig)%orientation(:,2) * lybyd_small

					small_config(iconfig)%node(5,:) = small_config(iconfig)%node(1,:) + &
																	small_config(iconfig)%orientation(:,3) * lybyd_small
					small_config(iconfig)%node(6,:) = small_config(iconfig)%node(2,:) + &
																	small_config(iconfig)%orientation(:,3) * lybyd_small
					small_config(iconfig)%node(7,:) = small_config(iconfig)%node(3,:) + &
																	small_config(iconfig)%orientation(:,3) * lybyd_small
					small_config(iconfig)%node(8,:) = small_config(iconfig)%node(4,:) + &
																	small_config(iconfig)%orientation(:,3) * lybyd_small
				endif
			enddo
		endif

		filename2 = "small_confing_coord.dat"
		open(unit=1, file=trim(filename2), status="replace", action="write")
		do iconfig=1, nsmall_box
			write (1,*) "variables=V1,V2,V3"
			write (1,*) "zone"
			write (1,*) "I=2, J=2, K=2, F=point"
			write (1,"(3d15.7)") small_config(iconfig)%node(1,:)
			write (1,"(3d15.7)") small_config(iconfig)%node(2,:)
			write (1,"(3d15.7)") small_config(iconfig)%node(4,:)			
			write (1,"(3d15.7)") small_config(iconfig)%node(3,:)
			write (1,"(3d15.7)") small_config(iconfig)%node(5,:)
			write (1,"(3d15.7)") small_config(iconfig)%node(6,:)
			write (1,"(3d15.7)") small_config(iconfig)%node(8,:)			
			write (1,"(3d15.7)") small_config(iconfig)%node(7,:)		
		enddo
		close(1)

		write (*,*) "LENGTH SEGMENTS = ", lybyd_small
		write (*,*) "NUMBER OF SMALL BOXES = ", nsmall_box

		do iconfig=1, nsmall_box
			pcount = 0
			do ibody=1, nbody
				vec1(:) = part(nzone, ibody, :)
				vec2(:) = vec1(:) - small_config(iconfig)%center(:)
				do idim=1, ndim
					if (vec2(idim)>rmax) then
						vec2(idim) = vec2(idim)-2*rmax
					elseif (vec2(idim)<-rmax) then
						vec2(idim) = vec2(idim)+2*rmax
					endif
				enddo

				do idim=1, ndim
					vec3(:) = small_config(iconfig)%orientation(:,idim)
					if( abs(dot_product(vec2,vec3)-lybyd_small/2)/sqrt(dot_product(vec3,vec3)) <= lybyd_small .and. &
						 abs(dot_product(vec2,vec3)+lybyd_small/2)/sqrt(dot_product(vec3,vec3)) <= lybyd_small) then

						inbox=.true.
					else
						inbox = .false.
						exit
					endif
				enddo
				if (inbox) pcount = pcount + 1
			enddo

			small_config(iconfig)%nbody = pcount

			if (small_config(iconfig)%nbody>0) then
				allocate(small_config(iconfig)%part(pcount))

				pcount = 0
				do ibody=1, nbody
					vec1(:) = part(nzone, ibody, :)
					vec2(:) = vec1(:) - small_config(iconfig)%center(:)
					do idim=1, ndim
						if (vec2(idim)>rmax) then
							vec2(idim) = vec2(idim)-2*rmax
						elseif (vec2(idim)<-rmax) then
							vec2(idim) = vec2(idim)+2*rmax
						endif
					enddo

					do idim=1, ndim
						vec3(:) = small_config(iconfig)%orientation(:,idim)
						if( abs(dot_product(vec2,vec3)-lybyd_small/2)/sqrt(dot_product(vec3,vec3)) <= lybyd_small .and. &
							 abs(dot_product(vec2,vec3)+lybyd_small/2)/sqrt(dot_product(vec3,vec3)) <= lybyd_small) then

							inbox=.true.
						else
							inbox = .false.
							exit
						endif
					enddo

					if (inbox) then
						pcount = pcount + 1
						small_config(iconfig)%part(pcount)%rad = part(nzone, ibody,ndim+1)
						small_config(iconfig)%part(pcount)%xc(:) = vec2(:)

					endif
				enddo

				if (pcount/=small_config(iconfig)%nbody) then
					write (*,*) "MISMATCH IN NUMBER OF PARTICLES, ", small_config(iconfig)%nbody, pcount
				endif
			endif
		enddo

		write (*,*) "SUM OF PARTICLES = ", sum(small_config(:)%nbody), nbody


!		filename2 = "small_confing_coord.dat"
!		open(unit=1, file=trim(filename2), status="replace", action="write")

		rmax = lybyd_small/2
		do iconfig=1, nsmall_box
#if 0
			write (1,*) "variables=V1,V2,V3"
			write (1,*) "zone"
			write (1,*) "I=2, J=2, K=2, F=point"
			write (1,"(3d15.7)") small_config(iconfig)%node(1,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(2,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(4,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(3,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(5,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(6,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(8,:)-small_config(iconfig)%center(:)
			write (1,"(3d15.7)") small_config(iconfig)%node(7,:)-small_config(iconfig)%center(:)
#endif

			nbody = small_config(iconfig)%nbody
			if (rotated) then
				do dim1=1, ndim
					do dim2=1, ndim
						rotation_tensor(dim1,dim2) = small_config(iconfig)%orientation(dim2,dim1)
					enddo
				enddo

				do ibody=1, nbody
					vec1(:) = small_config(iconfig)%part(ibody)%xc(:)
					call rotate_vec_backward(rotation_tensor, vec1, vec2, ndim)
					small_config(iconfig)%part(ibody)%xc(:) = vec2(:)
				enddo
			endif

			do ibody=1, small_config(iconfig)%nbody
				small_config(iconfig)%part(ibody)%xc(:) = small_config(iconfig)%part(ibody)%xc(:) + lybyd_small/2
			enddo

			ibody = 1
			do !ibody=1, small_config(iconfig)%nbody-1
				jbody = ibody+1
				do !jbody=1, small_config(iconfig)%nbody
					do idim=1, ndim
						vec1(idim) = abs( small_config(iconfig)%part(ibody)%xc(idim)-small_config(iconfig)%part(jbody)%xc(idim) )
						if (vec1(idim)>rmax) vec1(idim) = 2*rmax - vec1(idim)
					enddo
					small_max_overlap = dia_phys-sqrt(dot_product(vec1, vec1))

					if (small_max_overlap>rad_factor*dia_phys) then
						write (*,*) "REMOVING A PARTICLE, OVERLAP = ", small_max_overlap
						do kbody=jbody+1, small_config(iconfig)%nbody
							small_config(iconfig)%part(kbody-1)%xc(:) = small_config(iconfig)%part(kbody)%xc(:)
							small_config(iconfig)%part(kbody-1)%rad = small_config(iconfig)%part(kbody)%rad
						enddo
						small_config(iconfig)%part( small_config(iconfig)%nbody )%xc(:) = zero
						small_config(iconfig)%part( small_config(iconfig)%nbody )%rad = zero

						small_config(iconfig)%nbody = small_config(iconfig)%nbody - 1
					else
						jbody = jbody+1
					endif
					if (jbody>small_config(iconfig)%nbody) then
						exit
					endif
				enddo
				ibody = ibody+1
				if (ibody>small_config(iconfig)%nbody) then
					exit
				endif
			enddo

			small_config(iconfig)%volfrac = small_config(iconfig)%nbody*pi/6/lybyd_small**3
			write (*,"(1a,1i,1a,1i,1f8.4)") "# OF PARTICLES, VOLUME FRACTION IN BOX", iconfig, " = ", small_config(iconfig)%nbody, small_config(iconfig)%volfrac
		enddo
!		close(1)

		nbins = nbins / (lybyd/lybyd_small)
		allocate(gofr_small(nsmall_box,nbins), gofr_small_avg(nbins), gofr_small_var(nbins))
		allocate(rad_bin(nbins))

		gofr_small = zero
		gofr_small_avg = zero
		gofr_small_var = zero


		filename1 = "small_configs.dat"
		filename2 = "small_gofr.dat"
		filename4 = "small_cluster"

		filename5 = "small_cluster_projection.dat"
		filename6 = "small_cluster_area.dat"

!		open (unit=1, file=trim(filename1), status="replace", action="write")
		open (unit=2, file=trim(filename2), status="replace", action="write")

		open (unit=500, file=trim(filename5), status="replace", action="write")
		open (unit=600, file=trim(filename6), status="replace", action="write")

		first_pass = .true.
		do iconfig=1, nsmall_box
			nbody = small_config(iconfig)%nbody

			allocate(config_tmp(nbody,ndim), rad_tmp(nbody))
			allocate(contact(nbody,nbody))

			contact = .false.

			do ibody=1, nbody
				config_tmp(ibody,:) = small_config(iconfig)%part(ibody)%xc(:)
				rad_tmp(ibody) = small_config(iconfig)%part(ibody)%rad
			enddo

			small_max_overlap = max_overlap
			call calculate_gofr_homog(nbody, config_tmp, contact, int(lybyd_small), int(lybyd_small+1), nbins, .true., gofr_small(iconfig,1:nbins), rad_bin, small_max_overlap)

!			if (small_max_overlap>max_overlap) write (*,"(1a,1d15.7)") "max overlap = ", small_max_overlap

			call find_clusters(config_tmp, contact, nbody, lybyd_small, filename4, .true.)

!			write (1,*) "zone"
!			do ibody=1, nbody
!				write (1,"(4d15.7,1i)") small_config(iconfig)%part(ibody)%xc(1:ndim), small_config(iconfig)%part(ibody)%rad, aggregate_num(ibody)
!			enddo

			call output_particle_config1(config_tmp, aggregate_num, rad_tmp, nbody, filename1)
			first_pass = .false.

			write (2,*) "zone"
			do i=0, nbins
!				if (rad_bin(i)<dia_phys) gofr_small(iconfig,i)=zero
				write (2,"(2d15.7)") rad_bin(i), gofr_small(iconfig,i)
			enddo

#if 0
			do icluster=1, ncluster
				if (cluster_part_min<=cluster(icluster)%num .and. cluster(icluster)%num<=cluster_part_max) then
					write (600,*) "zone"
					surface(idim) = zero
					cord(idim) = 0

					allocate(cluster2d(cluster(icluster)%num,2))
					cluster2d = zero

					do idim=1, ndim
						write (500,*) "zone"

						do i=1, cluster(icluster)%num
							if (idim==1) then
								dim(1) = 2
								dim(2) = 3
							elseif (idim==2) then
								dim(1) = 3
								dim(2) = 1
							elseif (idim==3) then
								dim(1) = 1
								dim(2) = 2
							endif

							do jdim=1, 2
								cluster2d(i,jdim) = config_tmp( cluster(icluster)%list(i),dim(jdim) )
							enddo
						enddo

						write (500,"(3d15.7)") ((cluster2d(i,:), rad_tmp(i)), i=1, cluster(icluster)%num)

						surface(idim) = pi*rad_tmp(1)**2 * cluster(icluster)%num

						do i=1, cluster(icluster)%num-1
							do j=i+1, cluster(icluster)%num
								do jdim=1, 2
									dist(jdim) = abs( cluster2d(j,jdim)-cluster2d(i,jdim) )

									if (jdim==1 .and. dist(jdim)>cord(idim)) cord(idim) = dist(jdim)

									if (dist(jdim)>rmax) dist(jdim) = 2*rmax-dist(jdim)
								enddo

								sep = sqrt( dot_product(dist,dist) )
									
						
								if (sep<2*rad_tmp(1)) then
									alpha = asin(sep/2 / rad_tmp(1))

									surface(idim) = surface(idim) - 0.5 * rad_tmp(1)**2 * (2*alpha-sin(2*alpha))
								endif
							enddo
						enddo
					enddo
					deallocate(cluster2d)

					write (600,"(12d15.7)") surface(:)/(pi*rad_tmp(1)**2), surface(1)/surface(2), surface(1)/surface(3), surface(2)/surface(3), &
										 & cord(:)/(2*rad_tmp(1)), cord(1)/cord(2), cord(1)/cord(3), cord(2)/cord(3)
				endif
			enddo
#endif

			deallocate(config_tmp, rad_tmp)
			deallocate(contact, aggregate_num)
		enddo
!		close(1)
		close(2)

		close(500)
		close(600)

!		i = iconfig-1
		call get_confin(nsmall_box, confint)

		call calc_avr_var(nsmall_box, small_config(1:nsmall_box)%volfrac , volfrac_avg, volfrac_var)
		volfrac_var = sqrt(volfrac_var / (nsmall_box*(nsmall_box-1)))
		write (*,"(1A,2f8.4)") "AVERAGE VOL_FRAC = ", volfrac_avg, volfrac_var

		gofr_var = zero
		gofr_2d  = zero
		do ibin=1, nbins
			call calc_avr_var(nsmall_box, gofr_small(1:nsmall_box,ibin), gofr_small_avg(ibin), gofr_small_var(ibin))

			gofr_small_var(ibin) = sqrt(gofr_small_var(ibin) / (nsmall_box*(nsmall_box-1)))
		enddo

		filename3 = "small_gofr_avg.dat"
		open (unit=3, file=trim(filename3), status="replace", action="write")
		write (3,"(3d15.7)") ((rad_bin(ibin), gofr_small_avg(ibin), gofr_small_var(ibin)), ibin=1,nbins)		
		close(3)

		deallocate(gofr_small, gofr_small_avg, gofr_small_var,rad_bin)
		deallocate(part)

		write (*,*) "LEAVING COHISIVE_PARTS..."
		call screen_separator(30,'-')
	end subroutine cohesive_particles

#if 0
	subroutine simulated_annealing
		implicit none

		integer :: i, j, k, m, ibin, i_accept_total, j_accept, i_not_accept
!		integer :: nx1, ny1, nz1
		character*50 filename1, filename2, filename3, filename4, filename5
		logical :: accept
		real(prcn) :: obj, obj_old, obj_var, object0, temperature0
		integer :: nvar, nvar1, nvar2
		real(prcn), allocatable :: obj_array(:)

		first_pass = .true.
		temperature = temp_init
		call screen_separator(30,'^')
		write (*,*) "IN SIMULATED ANNEALING"

		if (post_no_flow_mem_alloc) goto 100

		filename1 = trim(run_name)//"_part_confing_SA.dat"
		filename2 = trim(run_name)//"_gofr_2d.dat"
		filename3 = trim(run_name)//"_cofr_2d.dat"
		filename4 = trim(run_name)//"_cofrf_2d.dat"

!		allocate(xc_conf(nbody,ndim))
!		nx1 = 6
!		ny1 = 6
!		nz1 = 6
!		nbody = nx1*ny1*nz1

		allocate(color(nbody,ndim))


		call initialize_sa

		if (isa==1) then
			call initialize_lattice
		else
			xc_conf = xc
		endif

		call output_particle_config1(xc_conf, int(color(:,1)), rad, nbody, filename1)
!		call uni_dist(ran_num)
!		idim = 0
!		do k=1, nz1
!			do j=1, ny1
!				do i=1, nx1
!					idim = idim+1
!					rad(idim) = radbdy(1)
!
!					xc_conf(idim,1) = nx/2 -nx1/2*rad(idim) + (mod(i-1,nx1)+1) * rad(idim) !+ radbdy(idim) * ran_num(idim*2-1)
!					xc_conf(idim,2) = my/2 -ny1/2*rad(idim) + (mod(j-1,ny1)+1) * rad(idim) !+ radbdy(idim) * ran_num(idim*2)
!					xc_conf(idim,3) = mz/2 -nz1/2*rad(idim) + (mod(k-1,nz1)+1) * rad(idim) !+ radbdy(idim) * ran_num(idim*2)
!
!					xc_conf(idim,1) = ran_num(idim*ndim-2) * nx
!					xc_conf(idim,2) = ran_num(idim*ndim-1) * my
!					xc_conf(idim,3) = ran_num(idim*ndim)   * mz
!
!					color(idim,1) = 256. * i/nx1
!					color(idim,2) = 256. * j/ny1
!					color(idim,3) = 256. * k/nz1
!				enddo
!			enddo
!		enddo

		obj_old = large_number

		xc_conf_old = xc_conf

		temperature0 = temperature
		i_accept_total = 0
		i_not_accept = 0
		do i=1, sa_max_itr

!			if (mod(i-1,100)==0) 
!write (*,"(1a,2i,1d15.7)") "SA ITR, ACCEPTED, TEMP = ",i, i_accept_total, temperature


!			if (isa==1) then
!				sa_in_itr = nbody !*ndim
!			else 
!				sa_in_itr = 1
!			endif

			if (.not.allocated(obj_array)) allocate(obj_array(sa_in_itr))

			j_accept = 0
			obj_array = zero
			do j=1, sa_in_itr

				xc_conf_old = xc_conf
				if (isa==1) call refinement_single_part(xc_conf, rad, nbody)
				call dem(xc_conf, color, rad, nbody)
!				call calc_objective_function_single_part(obj)
!				call validate(temperature, obj, obj_old, accept)
				accept = .true.

				if (accept==.true.) then
					i_accept_total = i_accept_total+1
					j_accept = j_accept+1
					obj_array(j_accept) = obj

					if (i_accept_total==1) object0 = obj
				else
					xc_conf = xc_conf_old
				endif
			enddo

			if (j_accept>0) then
!				call calc_avr_var(j_accept, obj_array(1:j_accept), obj, obj_var)
!				obj_var = sqrt(obj_var)
!
!				obj_old = obj
!				call output_error
!
!				call output_function(gofr_2d, gofr_var, rad_bin, nbins+1, filename2)
!
				if (first_pass) first_pass = .false.
				call output_particle_config1(xc_conf, int(color(:,1)), rad, nbody, filename1)

!				call reduce_temperature
				i_not_accept = 0
			else
				i_not_accept = i_not_accept + 1
				if (mod(i_not_accept,skip_num)==0) then
					write (*,*) "NO FURTHER REDUCTION. EXIT SA..."
					exit
				endif

				if (mod(i_not_accept,20)==0) then
					rad_ratio = rad_ratio * rad_reduction_ratio
					write (*,*) "REDUCING RAD_RATIO TO = ", rad_ratio
				endif
			endif

			write (*,"(1i,4d15.7)") i, dble(j_accept)/sa_in_itr, obj/object0, obj_var/object0, temperature/temperature0
#if 0

!			if (obj<sa_error) then
!				write (*,*) "SA CONVERGED. ERROR = ", obj
!				exit
!			endif

			if (rad_ratio<sa_error) then
				write (*,*) "PARTICLE RELOCATION IS NEGLIGIBLE, i.e. ", rad_ratio/rad(1)
				write (*,*) "EXIT SA.."
				exit
			endif

			if (mod(i,50)==0) then
				first_pass = .true.
				filename5 = trim(run_name)//"_gofr_2d_avg.dat"
				call output_function(gofr_2d, gofr_var, rad_bin, nbins+1, filename5)
				first_pass = .false.
			endif
#endif
		enddo

		write (*,"(1a,2i)") "END OF SA. ITR & ACCEPTED ITR =", i, i_accept_total


		WRITE(*,*) 'WRITING IMMERSED BODY DATA'
		open(unit=1,file=TRIM(RUN_NAME)//'_sphr_center.inp',status='replace', Action='write')
		write (1,*) "zone"
!		do i = 1, nbody
!			if (xc_conf(i,2)/dbydx<=lybyd/4 .or. xc_conf(i,2)/dbydx>=lybyd * 3/4) write (1,"(4d15.7,1i)") xc_conf(i,:)/dbydx, rad(i)/dbydx, 1
!		enddo
!		close(1)

		do i=1, 1
			do j=1, 1
				do k=1, 1
					do m=1, nbody
						write (1,"(4d15.7,1i)") xc_conf(m,1)/dbydx+(i-1)*lybyd, &
					&									xc_conf(m,2)/dbydx+(j-1)*lybyd, &
					&									xc_conf(m,3)/dbydx+(k-1)*lybyd, rad(m), 1
					enddo
				enddo
			enddo
		enddo

		write (1,*) "zone"
		do i=1, 5
			do j=1, 2
				do k=1, 2
					do m=1, nbody
!						if (xc_conf(m,2)/dbydx + (j-1)*lybyd <= 3*lybyd/4 .or. xc_conf(m,2)/dbydx + (j-1)*lybyd >= 3*lybyd * 3/4) &
						write (1,"(4d15.7,1i)") xc_conf(m,1)/dbydx+(i-1)*lybyd, &
					&									xc_conf(m,2)/dbydx+(j-1)*lybyd, &
					&									xc_conf(m,3)/dbydx+(k-1)*lybyd, rad(m), 1
					enddo
				enddo
			enddo
		enddo
		close(1)

#if 0
		first_pass = .true.
		filename2 = trim(run_name)//"_gofr_2d_avg.dat"
		call output_function(gofr_2d, gofr_var, rad_bin, nbins+1, filename2)
#endif
		call finalize

100	continue
#if 0
		if (post_no_flow_mem_alloc) then

			filename2 = "_gofr_2d_avg.dat"
			filename3 = "_cofr_2d_avg.dat"
			filename4 = "_cofrf_2d_avg.dat"

			nvar  = 3
			nvar1 = 1
			nvar2 = 2
			i = nbins+1

			call mis_average(nvar, nvar1, nvar2, filename2, i)
			call mis_average(nvar, nvar1, nvar2, filename3, i)
			call mis_average(nvar, nvar1, nvar2, filename4, i)
		endif
#endif
	contains
		subroutine output_error
			implicit none

			character*50 filename1

			filename1 = trim(run_name)//"_temperature.dat"
			if (first_pass) then
				open (unit=1, file=trim(filename1), status="replace")
			else
				open (unit=1, file=trim(filename1), status="old", position="append")
			endif

			write (1,"(1i,4d15.7)") i, dble(j_accept)/sa_in_itr, obj_old/object0, obj_var/object0, temperature
			close(1)
		end subroutine output_error


		subroutine reduce_temperature
			implicit none

				temperature = temperature*cooling_rate

		end subroutine reduce_temperature

		subroutine validate(temp, obj, obj_old, accept)
			implicit none
			real(prcn), intent(in) :: temp, obj, obj_old
			logical, intent(out) :: accept

			real(prcn) :: tmp

			if (obj<obj_old) then
				accept = .true.
			else
				call rand_out(tmp)
				if (  exp(-(obj-obj_old)/temp)>tmp ) then
					accept = .true.
				else
					accept = .false.
				endif
			endif			
		end subroutine validate

		subroutine initialize_lattice
			implicit none

			integer :: nx1, ny1, nz1, nb_tmp
			real(prcn) :: dx_g, dy_g, dz_g, tmp
			integer :: i,j,k,idim

			nx1 = dble(nbody)**(dble(1)/3)
			ny1 = nx1
			nz1 = nx1

			nb_tmp = 0
			do k=1, nz1
				do j=1, ny1
					do i=1, nx1
						nb_tmp = nb_tmp+1
						rad(nb_tmp) = dbydx * half

						xc_conf(nb_tmp,1) = dble(i-1)/nx1*nx! * rad(nb_tmp)
						xc_conf(nb_tmp,2) = dble(j-1)/ny1*my! * rad(nb_tmp)
						xc_conf(nb_tmp,3) = dble(k-1)/nz1*mz! * rad(nb_tmp)

						color(nb_tmp,1) = 256. * dble(i)/nx1
						color(nb_tmp,2) = 256. * dble(j)/ny1
						color(nb_tmp,3) = 256. * dble(k)/nz1
					enddo
				enddo
			enddo

			if (nb_tmp<nbody) then

				do
					nb_tmp = nb_tmp + 1
					rad(nb_tmp) = dbydx * half

					call rand_out(tmp)
					xc_conf(nb_tmp,1) = tmp * nx
					color(nb_tmp,1) = 256. * tmp

					call rand_out(tmp)
					xc_conf(nb_tmp,2) = tmp * my
					color(nb_tmp,2) = 256. * tmp

					call rand_out(tmp)
					xc_conf(nb_tmp,3) = tmp * mz
					color(nb_tmp,3) = 256. * tmp

					if (nb_tmp==nbody) exit
				enddo
			endif

		end subroutine initialize_lattice
	end subroutine simulated_annealing

	subroutine initialize_arrays
		implicit none

		real(prcn) :: rmin, rmax, dr
		integer :: i, nplane

		if (allocated(rad_bin)) deallocate(rad_bin)
		allocate(rad_bin(0:nbins))
		allocate(gofr_2d(0:nbins), gofr_var(0:nbins))
		allocate(cofr_2d(0:nbins), cofr_var(0:nbins))
		allocate(gofr_3d(0:nbins), cofr_3d(0:nbins))
		allocate(cofrf_2d(0:nbins), cofrf_var(0:nbins))

		allocate(gofr_3d_main(0:nbins))
		allocate(cofr_3d_main(0:nbins))
		allocate(gofr_2d_main(0:nbins))
		allocate(cofr_2d_main(0:nbins))

		nplane = nx+my+mz
		allocate(cofr(0:nbins,nplane), cofrf(0:nbins,nplane))
		allocate(gofr(0:nbins,nplane))

		allocate(xc_conf(nbody,ndim), xc_conf_old(nbody,ndim), color(nbody,ndim))
		allocate(rad(nbody))
		allocate(center_changed(nbody))

		rad_bin = zero
		gofr_2d  = zero
		gofr_3d  = zero
		cofr_2d = zero
		cofr_3d = zero
		cofrf_2d = zero
		gofr = zero
		cofr = zero
		cofrf = zero

		gofr_3d_main = zero
		cofr_3d_main = zero
		gofr_2d_main = zero
		cofr_2d_main = zero

		xc_conf = zero
		color(nbody,ndim) = zero
		rad = zero
		center_changed = .false.


		rmin = zero
		rmax = dble(lybyd) * half
		dr = (rmax-rmin)/(nbins)
		do i=0, nbins
			rad_bin(i) = i*dr
		enddo
	end subroutine initialize_arrays

	subroutine initialize_sa
		implicit none

		real(prcn) :: rmin, rmax, dr
		integer :: i, ifile, iline, io
		logical :: filexist

		integer :: nplane, nline, nvar, maxline=1000
		real(prcn), allocatable :: rad_tmp(:), array_tmp(:,:), var_tmp(:)
		real(prcn) :: tmp, ratio
		character*50 filename1


		if (isa==1) then
			write (*,*) "RUNNING SIMULATED ANNEALING"

			nbody = phiavg*6/pi*lybyd**3
			write (*,*) "THE NUMBER OF PARTICLES ARE ", nbody
			write (*,*) "THE EXACT VOLUME FRACTION =  ", nbody*pi/6/lybyd**3

			nx = lybyd*dbydx
			my = nx
			mz = nx
		endif

		call initialize_arrays

!		nbins = (lybyd*dbydy/7.5)*100
		write (*,*) "NBINS = ", nbins

		do ifile=1, 4
			if (ifile==1) then
				filename1 = "NMIS_gofr_main.dat"
				nvar = 6
			elseif (ifile==2) then
				filename1 = "NMIS_cofr_main.dat"
				nvar = 6
			elseif (ifile==3) then
				filename1 = "NMIS_gofr_2d_main.dat"
				nvar = 6
			elseif (ifile==4) then
				filename1 = "NMIS_cofr_2d_main.dat"
				nvar = 6
			endif

			inquire(file=trim(filename1), exist=filexist)
			if (.not.filexist) then
				write (*,*) 'FILE "'//trim(filename1)//' DOES NOT EXIST, SKIPPING '//trim(filename1)
			else
				write (*,*) trim(filename1)//" IS AVAILABLE, ALLOCATING VARIABLES..."
				open (unit=1,file=trim(filename1),status="old",action="read")

				allocate(array_tmp(nvar,maxline))
				allocate(var_tmp(nvar))

				array_tmp = zero
				var_tmp   = zero

				iline = 0
				do
					read (1,*,iostat=io) var_tmp(:)

					if (io>0) then
						write (*,*) "CHECK INPUT. SOMETHING IS WRONG IN LINE ", iline+1
					elseif (io==0) then
						iline = iline + 1
						array_tmp(1:nvar, iline) = var_tmp(1:nvar)
					elseif (io<0) then
						write (*,*) "END OF FILE ", trim(filename1)
						write (*,*) "NUMBER OF INPUTS = ", iline
						exit
					endif
				enddo
				nline = iline
				close(1)

				do i=0, nbins
					do iline=1, nline-1
						if ( (array_tmp(1,iline)<=rad_bin(i)).and.(rad_bin(i)<=array_tmp(1,iline+1)) ) then
							ratio = (rad_bin(i)-array_tmp(1,iline)) / (array_tmp(1,iline+1)-array_tmp(1,iline))

							tmp = (one-ratio)*array_tmp(2,iline) + ratio*array_tmp(2,iline+1)
							if (ifile==1) then
								gofr_3d_main(i) = tmp
							elseif (ifile==2) then
								cofr_3d_main(i) = tmp
							elseif (ifile==3) then
								gofr_2d_main(i) = tmp
							elseif (ifile==4) then
								cofr_2d_main(i) = tmp
							endif
						endif
					enddo
				enddo

				deallocate(array_tmp, var_tmp)
			endif
		enddo

		rad_bin = rad_bin*dbydx
	end subroutine initialize_sa

	subroutine finalize
		implicit none
		character*50 filename

		deallocate(rad_bin)
		deallocate(gofr_2d, gofr_var, gofr_3d)
		deallocate(cofr_2d, cofr_var, cofr_3d)
		deallocate(cofrf_2d, cofrf_var)

		if (allocated(gofr_2d_main)) deallocate(gofr_2d_main)
		if (allocated(cofr_2d_main)) deallocate(cofr_2d_main)
		if (allocated(gofr_3d_main)) deallocate(gofr_3d_main)
		if (allocated(cofr_3d_main)) deallocate(cofr_3d_main)

!		if (isa==1) then
			FILENAME = TRIM(RUN_NAME)//'_sphere_config.rst'
			open(unit=1,file=FILENAME,form="unformatted",status="unknown")

			write (1) nbody
			write (1) 1
			write (1) nbody
			write (1) rad(1)*two
			write (1) lybyd, lybyd, lybyd

			write(1) one
			write(1) xc_conf(1:nbody, 1:3)
			write(1) rad(1:nbody)
			write(1) color(1:nbody, 1:3)
			close(1)
!		endif

		deallocate(xc_conf, xc_conf_old, color, rad)
	end subroutine finalize


	subroutine output_function(func, func_var, rad, n, name)
		implicit none
		integer, intent(in) :: n
		real(prcn), intent(in) :: func(n), func_var(n), rad(n)
		character*50, intent(in) :: name
		integer :: i

		if (first_pass) then
			open (unit=1, file=trim(name), status="replace")
		else
			open (unit=1, file=trim(name), status="old", position="append")
		endif
		write (1,*) "zone"
		write (1,"(1f8.4,2d15.7)") ((rad(i)/(dbydx), func(i), func_var(i)), i=1,n)
		close(1)
	end subroutine output_function

	subroutine refinement_single_part(xc, rad, nbody)
		implicit none

		integer, intent(in) :: nbody
		real(prcn), intent(in) :: rad(nbody)
		real(prcn), intent(inout) :: xc(nbody,ndim)
		
		integer :: j, ibody
		real(prcn) :: tmp

!		call rand_out(tmp)
!		ibody = tmp*nbody+1
!		if (ibody>nbody) ibody = nbody

do ibody=1, nbody
		do j=1, ndim
			call rand_out(tmp)
			tmp = tmp-half

			xc(ibody,j) = xc(ibody,j) + two*rad(ibody)*rad_ratio*tmp
		enddo
enddo
	end subroutine refinement_single_part

	subroutine dem(xc, color, rad, nbody)
		implicit none

		integer, intent(in) :: nbody
		real(prcn), intent(in) :: color(nbody,ndim), rad(nbody)
		real(prcn), intent(inout) :: xc(nbody,ndim)
		

		logical, allocatable :: contact_list(:,:)
		real(prcn), allocatable :: sf_tot(:,:), normal(:,:,:)

		integer :: i, j, k, count, idim
		integer(8), save :: dem_itr=0
		real(prcn) :: vec(ndim), rmax(ndim), l(ndim), dist1, dist2, minforce
		real(prcn) :: sconst=1
		logical :: converge
		character*50 filename1, filename2

		logical, save :: first_time=.true.


		if (first_time) then
			center_changed = .true.
			first_time = .false.
		else
			center_changed = .false.
		endif

		filename1 = "_part_config_dem.dat"

		filename2 = trim(run_name)//"_dem_error.dat"
		open (unit=2, file=trim(filename2), status="replace")
		write (2,*) "zone"

		allocate(contact_list(nbody,nbody))
		allocate(normal(nbody, nbody, ndim+1))
		allocate(sf_tot(nbody,ndim+1))

		l(1) = dble(nx)
		l(2) = dble(my)
		l(3) = dble(mz)

		rmax(1) = l(1)/2
		rmax(2) = l(2)/2
		rmax(3) = l(3)/2

		converge = .false.
		dem_itr = 0
		do
			dem_itr = dem_itr+1
			contact_list = .false.
			normal = zero
			do i=1, nbody
				do j=1, i
					if (i/=j) then
90						continue
						dist1 = zero
						do idim=1, ndim
							vec(idim) = xc(j,idim)-xc(i,idim)
							if (abs(vec(idim))>rmax(idim)) then
								if (vec(idim)>0) then
									vec(idim) = vec(idim)-l(idim)
								else
									vec(idim) = l(idim)+vec(idim)
								endif
							endif

							dist1 = dist1 + vec(idim)**2
						enddo
						dist1 = sqrt(dist1)

						if (dist1<post_small) then
							write (*,*) "TWO NODES AT THE SAME PLACE"
							xc(i,1) = xc(i,1)+1
							if (xc(i,1)>l(1)) xc(i,1) = xc(i,1) - l(1)
							goto 90
						endif

						dist2 = rad(i)+rad(j)

						if (dist2-dist1>post_small) then
							contact_list(i,j) = .true.
							normal(i,j,1:ndim) = -vec(:)/dist1
							normal(i,j,ndim+1) = dist2-dist1

							contact_list(j,i) = contact_list(i,j)
							normal(j,i,1:ndim) = -normal(i,j,1:ndim)
							normal(j,i,ndim+1) = normal(i,j,ndim+1)
						endif
					endif
				enddo
			enddo

!			if (mod(dem_itr-1,20)==0) write (*,"(1a,1i,1a,1d15.7)") "max of normals in itr", dem_itr, "=", maxval(normal(:,:,ndim+1))/radbdy(1)
			write (2,"(1i10,1d15.7)") dem_itr, sum(normal(:,:,ndim+1))/rad(1)

			if (maxval(normal(:,:,ndim+1))/rad(1)<post_error) converge = .true.
			if (converge) exit

			vec = zero
			sf_tot = zero
			do i=1, nbody
				do j=1,nbody
					if (contact_list(i,j)==.true.) then
							sf_tot(i,1:ndim) = sf_tot(i,1:ndim) + normal(i,j,1: ndim) * normal(i,j,ndim+1)/2 * sconst
					endif
				enddo
				sf_tot(i,ndim+1) = sqrt(sum(sf_tot(i,1:ndim)**2))

				if (sf_tot(i,ndim+1)>post_small) sf_tot(i,1:ndim) = sf_tot(i,1:ndim)/sf_tot(i,ndim+1)
			enddo
#if 0
			minforce = -one
			do i=1, nbody
				if (sf_tot(i,ndim+1)>post_small) then
					if (minforce<zero) then
						minforce = max(sf_tot(i,ndim+1), minforce)
					else
						minforce = min(sf_tot(i,ndim+1), minforce)
					endif
				endif
			enddo

			if (minforce>zero) then
				do i=1, nbody
					if (sf_tot(i,ndim+1)>post_small) sf_tot(i,ndim+1) = minforce
				enddo
			endif
#endif

#if 0
			write (*,*) "CHECKING BINARY FORCES. BEFORE:"
			do idim=1, ndim
				write (*,*) sum(sf_tot(:,idim)*sf_tot(:,4))
			enddo

			do i=1, nbody
				do j=1, nbody
					if (contact_list(i,j)) then
						if (abs(sf_tot(i,ndim+1)*dot_product(sf_tot(i,1:ndim),normal(i,j,1:ndim)))>normal(i,j,ndim+1)/2) then
							sf_tot(i,ndim+1) = normal(i,j,ndim+1)/2 / abs(dot_product(sf_tot(i,1:ndim),normal(i,j,1:ndim)))
						endif
					endif
				enddo
			enddo

			write (*,*) "CHECKING BINARY FORCES. AFTER:"
			do idim=1, ndim
				write (*,*) sum(sf_tot(:,idim)*sf_tot(:,4))
			enddo
#endif

			do i=1, nbody
				if (sf_tot(i,ndim+1)>post_small) then
					do idim=1, ndim
						xc(i,idim) = xc(i,idim) + sf_tot(i,idim)*sf_tot(i,ndim+1)
						if (xc(i,idim)>l(idim)) then
							xc(i,idim) = xc(i,idim)-l(idim)
						elseif (xc(i,idim)<0) then
							xc(i,idim) = l(idim)+xc(i,idim)
						endif
					enddo
					center_changed(i) = .true.
				endif
			enddo

!		do i=1, nbody
!			if (sf_tot(i,ndim+1)>post_small) then
!				vel(i,1:ndim) = vel(i,1:ndim) * vel(i,ndim+1) + sf_tot(i,1:ndim) * sf_tot(i,ndim+1) * delt
!				vel(i,ndim+1) = sqrt(sum(vel(i,1:ndim)**2))
!
!				if (vel(i,ndim+1)>post_small) vel(i,1:ndim) = vel(i,1:ndim)/vel(i,ndim+1)
!			else
!				vel(i,:) = zero
!			endif
!		enddo
!
!		do i=1, nbody
!			if (vel(i,ndim+1)>post_small) then
!				xc(i,1:ndim) = xc(i,1:ndim) + vel(i,1:ndim)*delt
!			endif
!		enddo
		enddo
		close(2)

		deallocate(contact_list, sf_tot, normal)
	end subroutine dem



	subroutine projection_2d_single_part
		implicit none

		integer :: i, j, k, m, dim1, start, end1, end2, end3, iplane, jplane, nplane, ibin
		integer :: left1, left2, left, right1, right2, right

		logical, allocatable :: fluid_ind_2d(:,:), modified(:)
		real(prcn) :: dist
		real(prcn) :: xc_2d(nbody,2), rad_2d(nbody), confint, rmin, rmax, dr
		integer :: nbody_2d
		integer :: nvar, nvar1, nvar2
		character*50 filename1, filename2, filename3

!		call screen_separator(30,'^')
!		write (*,*) "IN PROJECTION_2D..."


		nplane = nx+my+mz

!			filename1 = trim(run_name)//"_sphr_projection.dat"
!			open (unit=1, file=trim(filename1), status="replace")
!			write (1,"(1A)") "zone"
!			write (1,"(3I6)") nx, my, mz


		if (.not.allocated(modified)) allocate(modified(nplane))
		modified = .false.

#if 0
		do i=1, nbody
			if (.not.center_changed(i)) goto 10

			do dim1=1, ndim
				left1 = xc_conf(i, dim1)-rad(i)
				left2 = xc_conf_old(i, dim1)-rad(i)
				left = min(left1,left2)
				if (dim1==1) then
					if (left<1) left = left+nx
				elseif (dim1==2) then
					if (left<1) left = left+my
				elseif (dim1==3) then
					if (left<1) left = left+mz
				endif

				right1 = xc_conf(i, dim1)+rad(i)
				right2 = xc_conf_old(i, dim1)+rad(i)

				right = max(right1,right2)
				if (dim1==1) then
					if (right>nx) right = right-nx
				elseif (dim1==2) then
					if (right>my) right = right-my
				elseif (dim1==3) then
					if (right>mz) right = right-mz
				endif


				if (dim1==1) then
					start=0
					end1 = nx
				elseif (dim1==2) then
					start = nx
					end1 = my
				elseif (dim1==3) then
					start = nx+my
					end1 = mz
				endif

				if (left<right) then
					do j=left, right
						modified(start+j) = .true.
					enddo
				elseif (right<left) then
					do j=1, right
						modified(start+j) = .true.
					enddo

					do j=left, end1
						modified(start+j) = .true.
					enddo
				endif
			enddo
10			continue
		enddo
#endif


		iplane = 0
		jplane = 0
		start = 1
		do dim1=1, ndim
			if (dim1==1) then
				end1 = nx

				end2 = my
				end3 = mz
			elseif (dim1==2) then
				end1 = my

				end2 = nx
				end3 = mz
			elseif (dim1==3) then
				end1 = mz

				end2 = nx
				end3 = my
			endif

			allocate(fluid_ind_2d(end2,end3))

			do i=start, end1
				iplane = iplane+1

!				if (.not.modified(iplane)) goto 20

				if (isa==0) then
					if (mod(iplane, 20)==0) write (*,"(2(1A,1I))") "@ IPLANE = ", iplane, ", JLPANE = ", jplane
				endif

				xc_2d = zero
				rad_2d = zero
				nbody_2d = 0

				call project_particles

				if (nbody_2d==0) then
					write (*,*) "NO PROJECTED PARTICLE AT PLANE", dim1, i
				else
					!WRITING PROJECTED SURFACES
!						write (1,"(1A)") "zone"
!						write (1,"(3d15.7)") ((xc_2d(m,:), rad_2d(m)), m=1m nbody_2d)
				endif

				call calc_gofr_2d(nbody_2d, nbins, end2, end3, xc_2d(1:nbody_2d,1:2), rad_2d(1:nbody_2d), rad_bin, gofr(0:nbins,iplane))

#if 0
if (isa==0) then
				if (mod(iplane, skip_num)==0) then
					jplane = jplane+1

					if (dim1==1) then
						fluid_ind_2d(:,:) = fluid_atijk(i,:,:)
					elseif (dim1==2) then
						fluid_ind_2d(:,:) = fluid_atijk(:,i,:)
					elseif (dim1==3) then
						fluid_ind_2d(:,:) = fluid_atijk(:,:,i)
					endif

					call calc_cofr_2d(end2, end3, nbins, fluid_ind_2d, rad_bin, cofr(0:nbins,jplane), cofrf(0:nbins,jplane))
				endif
endif
#endif

20				continue
			enddo

			deallocate(fluid_ind_2d)
		enddo
!			close(1)

		!^^^^ OUTPUT FOR G(R) ^^^^^^
!			filename2 = trim(run_name)//"_gofr_2d_planes.dat"
!			open (unit=2, file=trim(filename2), status="replace")
!			do iplane=1, nplane
!				write (2,"(1A)") "zone"
!				do ibin=0, nbins
!!					if (gofr_2d(ibin,iplane)>post_small)
!					write (2,"(2d15.7)") rad_bin(ibin)/(radbdy(1)*2), gofr(ibin,iplane)
!				enddo
!			enddo
!			close(2)

!		write (*,*) "NUMBER OF PLANES USED TO FORM G(r) ENSEMBLE=", nplane

!		filename3 = trim(run_name)//"_gofr_2d_avg.dat"
!		if (first_pass) then
!			open (unit=3, file=trim(filename3), status="replace")
!		else
!			open (unit=3, file=trim(filename3), status="old", position="append")
!		endif
!		write (3,*) "zone"

		call get_confin(nplane, confint)
		gofr_var = zero
		do ibin=0, nbins
			call calc_avr_var(nplane, gofr(ibin,1:nplane), gofr_2d(ibin), gofr_var(ibin))

			gofr_var(ibin) = sqrt(gofr_var(ibin) / (nplane*(nplane-1)))
!				write (3,"(1f8.4,2d15.7)") rad_bin(ibin)/(radbdy(1)*2), gofr_2d(ibin), gofr_var(ibin)
		enddo
		close(3)

!		deallocate(gofr)
		!--------------------------

		!^^^^^ OUTPUT FOR C(R) ^^^^
!			filename2 = trim(run_name)//"_cofr_2d.dat"
!			open (unit=2, file=trim(filename2), status="replace")
!
!			do iplane=1, nplane
!				write (2,"(1A)") "zone"
!				do ibin=0, nbins
!					write (2,"(2d15.7)") rad_bin(ibin)/(radbdy(1)*2), cofr_2d(ibin,iplane)
!				enddo
!			enddo
!			close(2)
!		write (*,*) "NUMBER OF PLANES USED TO FORM C(r) ENSEMBLE=", jplane

!			filename3 = trim(run_name)//"_cofr_2d_avg.dat"
!			if (first_pass) then
!				open (unit=3, file=trim(filename3), status="replace")
!			else
!				open (unit=3, file=trim(filename3), status="old", position="append")
!			endif
!			write (3,*) "zone"

#if 0
if (isa==0) then
		call get_confin(jplane, confint)

		cofr_var = zero
		cofrf_var = zero

!		write (*,*) "COMPUTING AVERAGE OF TWO-POINT CORRELATION OF THE INDICATOR FIELD..."

		do ibin=0, nbins
			call calc_avr_var(jplane, cofr(ibin,1:jplane), cofr_2d(ibin), cofr_var(ibin))
			call calc_avr_var(jplane, cofrf(ibin,1:jplane), cofrf_2d(ibin), cofrf_var(ibin))

			cofr_var(ibin) = sqrt(cofr_var(ibin) / (jplane*(jplane-1)))
			cofrf_var(ibin) = sqrt(cofrf_var(ibin) / (jplane*(jplane-1)))

!				if (cofr_2d_avg(ibin)>post_small) then
!				write (3,"(1f8.4,2d15.7)") rad_bin(ibin)/(radbdy(1)*2), cofr_2d(ibin), cofr_var(ibin)
!				endif
		enddo
!		close(3)
		!-------------------------
!		deallocate(cofr, cofrf)
endif
#endif
!		call screen_separator(30,'-')

	contains

		subroutine project_particles
			implicit none
			!FINDING THE PROJECTION OF PARTICLES ON EACH PLANE				
			do m=1, nbody
				if (dim1==1) then
					dist = min(abs(xc_conf(m,1)-i), abs(xc_conf(m,1)-nx-i), abs(xc_conf(m,1)+nx-i))
				elseif (dim1==2) then
					dist = min(abs(xc_conf(m,2)-i), abs(xc_conf(m,2)-my-i), abs(xc_conf(m,2)+my-i))
				elseif (dim1==3) then
					dist = min(abs(xc_conf(m,3)-i), abs(xc_conf(m,3)-mz-i), abs(xc_conf(m,3)+mz-i))
				endif

				if (dist<rad(m)) then
					nbody_2d = nbody_2d+1
					if (dim1==1) then
						xc_2d(nbody_2d,1) = xc_conf(m,2)
						xc_2d(nbody_2d,2) = xc_conf(m,3)
					elseif (dim1==2) then
						xc_2d(nbody_2d,1) = xc_conf(m,1)
						xc_2d(nbody_2d,2) = xc_conf(m,3)
					elseif (dim1==3) then
						xc_2d(nbody_2d,1) = xc_conf(m,1)
						xc_2d(nbody_2d,2) = xc_conf(m,2)
					endif
					rad_2d(nbody_2d) = sqrt(rad(m)**2-dist**2)
				endif
			enddo
		end subroutine project_particles
	end subroutine projection_2d_single_part


	subroutine calc_objective_function_single_part(objective_value)
		implicit none

		real(prcn), intent(out) :: objective_value

		real(prcn) :: l2_gofr2d, l2_cofr2d, l2_gofr3d, l2_cofr3d
		character*50 filename1

		call projection_2d_single_part

		objective_value = zero

		if (allocated(gofr_2d_main)) call l2_norm(gofr_2d, gofr_2d_main, rad_bin, nbins, l2_gofr2d)
		objective_value = objective_value + objective_coef1*l2_gofr2d

!		if (allocated(cofr_2d_main)) call l2_norm(cofr_2d, cofr_2d_main, rad_bin, nbins, l2_cofr2d)
!		objective_value = objective_value + objective_coef2*l2_cofr2d

!		if (allocated(gofr_3d)) call l2_norm(gofr_3d, gofr_3d_main, rad_bin, nbins, l2_gofr3d)
!		if (allocated(cofr_3d)) call l2_norm(cofr_3d, cofr_3d_main, rad_bin, nbins, l2_cofr3d)		
	contains
		subroutine l2_norm(main, temp, bin, nbin, norm)
			implicit none
			integer, intent(in) :: nbin
			real(prcn), intent(in) :: main(0:nbin), temp(0:nbin), bin(0:nbin)
			real(prcn), intent(out) :: norm

			integer :: i

			norm = zero

			do i=0, nbin
!				if (bin(i)<func_rad_ratio*rad(1)) then
!				if (main(i)>one) then		
					if (main(i)>post_small.and.temp(i)>post_small) then
						norm = norm + (main(i)-temp(i))**2
					endif
!				endif
			enddo
			norm = sqrt(norm)
		end subroutine l2_norm
	end subroutine calc_objective_function_single_part



	subroutine indicator_output
		implicit none

		integer :: i,j,k
		character*50 filename1

		if (.not.post_no_flow_mem_alloc) then
			filename1 = trim(run_name)//"_indicator.dat"
			open (unit=1, file=trim(filename1), status="replace")
			write (1,"(3I)") nx, my, mz

			do k=1, mz
				do j=1, my
					do i=1, nx
						if (fluid_atijk(i,j,k)) then
							write (1,"(1I2)") 1
						else
							write (1,"(1I2)") 0
						endif
					enddo
				enddo
			enddo
			close (1)
		endif
	end subroutine indicator_output



	subroutine calc_gofr_2d(nbody, nbin, nx, ny, xc_2d, rad_2d, rad, gr)
		implicit none
		integer, intent(in) :: nbody, nbin, nx, ny
		real(prcn), intent(in) :: xc_2d(nbody,2), rad_2d(nbody), rad(0:nbin)
		real(prcn), intent(out) :: gr(0:nbin)

		integer :: i, j, ibin
		real(prcn) :: dist, dist1, dist2, dr, rmin, rmax, volbin
		integer(8) :: bin_sample(0:nbin)

		dr = rad(1)-rad(0)
		rmin = rad(0)
		rmax = rad(nbin)
		bin_sample = 0
		gr = zero

		do i=1, nbody-1
			do j=i+1, nbody
				dist1 = abs(xc_2d(i,1)-xc_2d(j,1))
				if (dist1>rmax) dist1 = 2*rmax-dist1

				dist2 = abs(xc_2d(i,2)-xc_2d(j,2))
				if (dist2>rmax) dist2 = 2*rmax-dist2

				dist = sqrt(dist1**2 + dist2**2)
				if (dist>rmax) goto 100

				ibin = dist/dr

				if (ibin>nbin) ibin = nbin

				gr(ibin) = gr(ibin)+1
				bin_sample(ibin) = bin_sample(ibin)+1
100			continue
			enddo
		enddo

		do ibin=0, nbin
			if (bin_sample(ibin)>0) then
				volbin = 2*pi*rad(ibin)*dr / (nx*ny)
				if (volbin>post_small) gr(ibin) = gr(ibin)/nbody / volbin *2/nbody
			endif
		enddo
	end subroutine calc_gofr_2d

	subroutine calc_cofr_2d(nx, ny, bin_num, index, rad, indicator, indicatorf)
		implicit none
		integer, intent(in) :: nx, ny, bin_num
		logical, intent(in) :: index(nx,ny)
		real(prcn), intent(in) :: rad(0:bin_num)
		real(prcn), intent(out) :: indicator(0:bin_num), indicatorf(0:bin_num)

		integer(8) :: bin_sample(0:bin_num), ibin

		type :: node_type
			real(prcn) :: x, y
			logical :: ind
		end type node_type

		type(node_type), allocatable :: node(:)

		real(prcn) :: dr, rmin, rmax, dist, dist1, dist2
		integer :: i, j, k, node_num

		indicator = zero
		indicatorf = zero
		bin_sample = zero

		node_num = nx*ny
		allocate(node(node_num))

		k = 0
		do i=1, nx
			do j=1, ny
				k = k+1
				node(k)%x = i
				node(k)%y = j
				node(k)%ind = index(i,j)
			enddo
		enddo

		if (k>node_num) then
			write (*,*) "NUMBER OF POINTERS DOES NOT MATCH THE NUMBER OF NODES."
			stop
		endif

		dr = rad(1)-rad(0)
		rmin = rad(0)
		rmax = rad(bin_num)

		do i=1, node_num-1
			do j=i, node_num
!				if (node(i)%ind.or.node(j)%ind) goto 100
				dist1 = abs(node(i)%x-node(j)%x)
				if (dist1>rmax) dist1 = 2*rmax-dist1

				dist2 = abs(node(i)%y-node(j)%y)
				if (dist2>rmax) dist2 = 2*rmax-dist2

				dist = sqrt(dist1**2 + dist2**2)
				if (dist>rmax) goto 100

				ibin = dist/dr+1

				if (ibin>bin_num) ibin = bin_num
				if (i==j) ibin = 0

				if ((node(i)%ind==.false.).and.(node(j)%ind==.false.)) then
					indicator(ibin) = indicator(ibin)+1
				elseif ((node(i)%ind==.true.).and.(node(j)%ind==.true.)) then
					indicatorf(ibin) = indicatorf(ibin)+1
				endif
				bin_sample(ibin) = bin_sample(ibin)+1
100			continue
			enddo
		enddo

		do i=0, bin_num
			if (bin_sample(i)>0) indicator(i) = indicator(i)/bin_sample(i)
			if (bin_sample(i)>0) indicatorf(i) = indicatorf(i)/bin_sample(i)
		enddo

		do i=bin_num, 0, -1
			if (indicator(0)>post_small) indicator(i) = indicator(i)/indicator(0)
			if (indicator(0)>post_small) indicatorf(i) = indicatorf(i)/indicatorf(0)
		enddo

		deallocate(node)
	end subroutine calc_cofr_2d

	subroutine calc_cofr_3d
		implicit none

		real(prcn) :: rmin, rmax, dr
		integer :: nvar, nvar1, nvar2
		character*50 filename1

		integer(8) :: bin_sample(0:nbins)

		type :: node_type
			real(prcn) :: x, y, z
			logical :: ind
		end type node_type

		type(node_type), allocatable :: node(:)

		real(prcn) :: dist, dist1, dist2, dist3
		integer :: ibin, i, j, k, l, node_num, inode


		call screen_separator(30,'^')
		write (*,*) "IN cofr_3d..."
		if (.not.post_no_flow_mem_alloc) then
			write (*,*) "TOTAL NUMBER OF NODES = ", nx*my*mz

			if (.not.allocated(rad_bin)) allocate(rad_bin(0:nbins))
			if (.not.allocated(cofr_3d)) allocate(cofr_3d(0:nbins))
			rad_bin  = zero
			cofr_3d = zero
			bin_sample = zero

			rmin = zero
			rmax = dble(max(nx,my,mz)) / 2
			dr = (rmax-rmin)/(nbins)
			do i=0, nbins
				rad_bin(i) = i*dr
			enddo

			node_num = nx*my*mz
			allocate(node(node_num))

			l = 0
			do i=1, nx
				do j=1, my
					do k=1, mz
						l = l+1
						node(l)%x = i
						node(l)%y = j
						node(l)%z = k
						node(l)%ind = fluid_atijk(i,j,k)
					enddo
				enddo
			enddo

			if (k>node_num) then
				write (*,*) "NUMBER OF POINTERS DOES NOT MATCH THE NUMBER OF NODES."
				stop
			endif


			inode = 0
			do i=1, node_num-1, skip_num
				inode = inode + 1
				if (mod(inode,10)==0) write (*,*) "NODE_NUM = ",inode
				do j=i, node_num

					dist1 = abs(node(i)%x-node(j)%x)
					if (dist1>rmax) dist1 = 2*rmax-dist1

					dist2 = abs(node(i)%y-node(j)%y)
					if (dist2>rmax) dist2 = 2*rmax-dist2

					dist3 = abs(node(i)%z-node(j)%z)
					if (dist3>rmax) dist3 = 2*rmax-dist3

					dist = sqrt(dist1**2 + dist2**2 + dist3**2)
					if (dist>rmax) goto 100

					ibin = dist/dr+1

					if (ibin>nbins) ibin = nbins
					if (i==j) ibin = 0

					if ((node(i)%ind==.false.).and.(node(j)%ind==.false.)) then
						cofr_3d(ibin) = cofr_3d(ibin)+1
					endif
					bin_sample(ibin) = bin_sample(ibin)+1
100					continue
				enddo
				if (mod(inode,1000)==0) then
					call output
				endif
			enddo
			deallocate(node)
			call output

			do i=0, nbins
				if (bin_sample(i)>0) cofr_3d(i) = cofr_3d(i)/bin_sample(i)
			enddo

			do i=nbins, 0, -1
				if (cofr_3d(0)>post_small) cofr_3d(i) = cofr_3d(i)/cofr_3d(0)
			enddo
		else
			nvar  = 2
			nvar1 = 1
			nvar2 = 1

			filename1 = "_cofr_3d.dat"
			i = nbins+1
			call mis_average(nvar, nvar1, nvar2, filename1, i)
		endif
		call screen_separator(30,'-')
	contains
		subroutine output
			implicit none

			integer :: i
			character*50 filename1

			real(prcn), allocatable :: ind_tmp(:)

			write (*,*) "GENERATING OUTPUT..."
			allocate(ind_tmp(0:nbins))

			do i=0, nbins
				if (bin_sample(i)>0) ind_tmp(i) = cofr_3d(i)/bin_sample(i)
			enddo

			do i=nbins, 0, -1
				if (ind_tmp(0)>1e-6) ind_tmp(i) = ind_tmp(i)/ind_tmp(0)
			enddo

			filename1 = trim(run_name)//"_cofr_3d.dat"
			open (unit=1, file=trim(filename1), status="replace")
!			write (1,"(1A)") "zone"
			do i=0, nbins
!				if (cofr_3d_tmp(ibin)>post_small) then
					write (1,"(2d15.7)") rad(i)/(dbydx), ind_tmp(i)
!				endif
			enddo
			close(1)

			deallocate(ind_tmp)
		end subroutine output
	end subroutine calc_cofr_3d

#endif

end module geometry_2d_3d
