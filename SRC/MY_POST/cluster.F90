module clustermod
	use precision  
	use constants  
	use global_data		
	use general_funcs

	implicit none

	type :: cluster_type
		integer :: num
		real(prcn) :: center(ndim), rg, rg_tensor(ndim,ndim), eigen_val(ndim), eigen_vec(ndim,ndim), force(ndim)

		integer, allocatable :: list(:)
	end type cluster_type
	type(cluster_type), allocatable :: cluster(:)

	integer, allocatable :: aggregate_num(:)
	integer :: ncluster

contains

	subroutine find_clusters(part, contact, nbody, lbyd, filename, make_output)
		implicit none
		integer, intent(in) :: nbody
		logical, intent(in) :: make_output
		real(prcn), intent(in) :: part(nbody,ndim)
		logical, intent(in) :: contact(nbody,nbody)
		real(prcn), intent(in) :: lbyd
		character*50, intent(in) :: filename
		integer :: i, j, icluster, length, index
		integer :: ifound, jfound
		real(prcn) :: rmax

		rmax = lbyd/2

		ncluster = 0
		if (allocated(cluster)) then
!				length = size(cluster,1)
			length = size(cluster)

!				if (ncluster/=length) write (*,*) "THE NUMBER OF CLUSTERS DOES NOT MATCH THE NCLUSTER!!!"

			do icluster=1, length
				if (allocated(cluster(icluster)%list)) deallocate(cluster(icluster)%list)
			enddo
			deallocate(cluster)
		endif
		ncluster = 0

		allocate(cluster(nbody))
		do i=1, nbody
			cluster(i)%num = 0
		enddo

		do i=1, nbody-1
			do j=i+1, nbody
				if (.not.contact(i,j)) goto 10
			
				ifound = 0
				jfound = 0
				do icluster=1, ncluster
					do index=1, cluster(icluster)%num
						if (cluster(icluster)%list(index)==i) then
							ifound = icluster
						endif

						if (cluster(icluster)%list(index)==j) then
							jfound = icluster
						endif
					enddo
				enddo

				if ((ifound>0.and.jfound>0) .and. (ifound==jfound) ) then

				elseif ((ifound>0.and.jfound>0) .and. (ifound/=jfound) ) then

					call merge_cluster(ifound, jfound)
					call del_cluster(jfound)

				elseif (ifound>0.and.jfound==0) then

					call add_particle(ifound, j)

				elseif (jfound>0.and.ifound==0) then

					call add_particle(jfound, i)

				else

					call add_cluster(i,j)

				endif
10				continue
			enddo
		enddo

!		write (*,*) "NCLUSTER = ", ncluster
		call count_particles
		call compute_rg

		if (allocated(u)) call compute_force_per_cluster

		if (make_output) then
			call cluster_distribution(filename)

			if (first_pass) then
				open (unit=10, file=trim(filename)//".dat", status="replace", action="write")
			else
				open (unit=10, file=trim(filename)//".dat", status="old", action="write", position="append")
			endif

			write (10,*) "zone"
			do icluster=1, ncluster
				write (10,"(1i,17d15.7)") cluster(icluster)%num, dble(cluster(icluster)%num)/nbody, cluster(icluster)%center(:)/dia_phys, &
							& cluster(icluster)%rg/(dia_phys/2), cluster(icluster)%rg_tensor(:,:) / cluster(icluster)%rg**2, &
							& cluster(icluster)%eigen_val(:)/cluster(icluster)%rg**2
			enddo
			close(10)
		endif

	contains

		subroutine add_cluster(i,j)
			implicit none

			integer, intent(in) :: i,j

			ncluster = ncluster+1
			allocate(cluster(ncluster)%list(10))
			cluster(ncluster)%num = 2
			cluster(ncluster)%list(1) = i
			cluster(ncluster)%list(2) = j
			cluster(ncluster)%list(3:10) = 0
		end subroutine add_cluster


		subroutine add_particle(iclus, i)
			implicit none

			integer, intent(in) :: iclus, i
			integer :: length
			integer, allocatable :: list_tmp(:)

			length = size(cluster(iclus)%list)

!			if (length==cluster(iclus)%num) then
!				write (*,*) "INCREASING THE CLUSTER LIST FOR CLUSTER", iclus
!				write (*,*) "THE NUMBER OF THE CURRENT LIST IS ", cluster(iclus)%num
!			endif

!			if (mod(cluster(iclus)%num,10)==0) then
			if (cluster(iclus)%num==length) then
				if (minval(cluster(iclus)%list(:))==0) then
					write (*,*) "INCORRECT ENTRY IN NEIGHBOR_LIST"
					write (*,*) "ICLUSTER & NUM = ", iclus, cluster(iclus)%num
					write (*,*) cluster(iclus)%list(:)
					read (*,*)
				endif

!				length = size(cluster(iclus)%list)

				allocate(list_tmp(length))
				list_tmp(:) = cluster(iclus)%list(:)

				deallocate(cluster(iclus)%list)
				allocate(cluster(iclus)%list(length+10))

				cluster(iclus)%list(:) = 0
				cluster(iclus)%list(1:length) = list_tmp(1:length)

				deallocate(list_tmp)
			endif

			cluster(iclus)%num = cluster(iclus)%num+1
			cluster(iclus)%list( cluster(iclus)%num ) = i
		end subroutine add_particle


		subroutine merge_cluster(iclust, jclust)
			implicit none

			integer, intent(in) :: iclust, jclust
			integer :: i

			do i=1, cluster(jclust)%num
				call add_particle(iclust, cluster(jclust)%list(i))
			enddo
		end subroutine merge_cluster


		subroutine del_cluster(iclus)
			implicit none
			integer, intent(in) :: iclus

			integer :: i, isize, imsize, inum

			do i=iclus+1, ncluster
				inum   = cluster(i)%num
				isize  = size(cluster(i)%list)
				imsize = size(cluster(i-1)%list)

				if (imsize<isize) then
					deallocate(cluster(i-1)%list)
					allocate(cluster(i-1)%list(isize))
				endif
				cluster(i-1)%list(:) = 0

				cluster(i-1)%num = inum
				cluster(i-1)%list(1:inum) = cluster(i)%list(1:inum)
			enddo

			cluster(ncluster)%num = 0
			deallocate(cluster(ncluster)%list)
			ncluster = ncluster-1
		end subroutine del_cluster

		subroutine count_particles
			implicit none

			integer :: i, j

			if (allocated(aggregate_num)) deallocate(aggregate_num)
			allocate(aggregate_num(nbody))
			aggregate_num = 1

			do i=1, ncluster
				do j=1, cluster(i)%num
					aggregate_num(cluster(i)%list(j)) = cluster(i)%num
				enddo
			enddo
		end subroutine count_particles

		subroutine compute_rg
			use general_funcs, only : eigen_val_vec
			implicit none

			integer :: i, j, dim1, dim2, idim
			real(prcn) :: vec(ndim)

			do i=1, ncluster
				cluster(i)%center = zero
				do j=1, cluster(i)%num
					if (j>1) then
						do idim=1, ndim
							vec(idim) = part( cluster(i)%list(j),idim ) - part( cluster(i)%list(1),idim )



							if (vec(idim)>rmax) then
								vec(idim) = vec(idim)-2*rmax
							elseif (vec(idim)<-rmax) then
								vec(idim) = 2*rmax + vec(idim)
							endif
						enddo
						vec(:) = vec(:) + part( cluster(i)%list(1),: )
					else
						vec(:) = part( cluster(i)%list(j),: )
					endif

					cluster(i)%center(:) = cluster(i)%center(:) + vec(:)
				enddo
				cluster(i)%center(:) = cluster(i)%center(:) / cluster(i)%num

				do idim=1,ndim
					if (cluster(i)%center(idim)>2*rmax) then
						cluster(i)%center(idim)=cluster(i)%center(idim)-2*rmax
					elseif (cluster(i)%center(idim)<0) then
						cluster(i)%center(idim)=cluster(i)%center(idim)+2*rmax
					endif
				enddo

				cluster(i)%rg = zero
				cluster(i)%rg_tensor = zero
				do j=1, cluster(i)%num
					do idim=1, ndim
						vec(idim) = cluster(i)%center(idim) - part( cluster(i)%list(j),idim )
						if (vec(idim)>rmax) then
							vec(idim) = vec(idim) - 2*rmax
						elseif (vec(idim)<-rmax) then
							vec(idim) = vec(idim) + 2*rmax
						endif
					enddo

					do dim1=1, ndim
						do dim2=1, ndim
							cluster(i)%rg_tensor(dim1,dim2) = cluster(i)%rg_tensor(dim1,dim2) + vec(dim1)*vec(dim2) / cluster(i)%num
						enddo
					enddo
!					cluster(i)%rg = cluster(i)%rg + dot_product(vec,vec)
				enddo
!				cluster(i)%rg = sqrt(cluster(i)%rg / cluster(i)%num)

				cluster(i)%rg = sqrt(cluster(i)%rg_tensor(1,1)+cluster(i)%rg_tensor(2,2)+cluster(i)%rg_tensor(3,3))
				call eigen_val_vec(cluster(i)%rg_tensor(:,:), ndim, cluster(i)%eigen_val(:), cluster(i)%eigen_vec(:,:))

				if (cluster(i)%num>2) then
					call cross_product(cluster(i)%eigen_vec(:,1),cluster(i)%eigen_vec(:,2), vec(:))
					if (dot_product(cluster(i)%eigen_vec(:,3), vec(:))<zero) cluster(i)%eigen_vec(:,3) = -cluster(i)%eigen_vec(:,3)
				endif
			enddo
		end subroutine compute_rg


		subroutine cluster_distribution(filename)
			implicit none
				
			character*50, intent(in) :: filename
			character*50 filename1
			integer :: nhist, i, cluster_part
			integer, allocatable :: hist(:)
			real(prcn), allocatable :: weight(:)

			nhist = maxval(cluster(:)%num)
			allocate(hist(nhist), weight(nhist))
			hist = 0
			weight = zero
			cluster_part = 0
			do i=1, ncluster
				j = cluster(i)%num
				cluster_part = cluster_part + j
				hist(j) = hist(j) + 1
				weight(j) = weight(j) + dble(cluster(i)%num)/nbody
			enddo			
			if (cluster_part<nbody) then
				hist(1) = nbody - cluster_part
				weight(1) = dble(hist(1))/nbody
			endif

			filename1 = trim(filename)//"_dist.dat"

			if (first_pass) then
				open (unit=11, file=trim(filename1), status="replace", action="write")
			else
				open (unit=11, file=trim(filename1), status="old", action="write", position="append")
			endif

			write (11,*) "zone"
			do i=1, nhist
				if (hist(i)>0) write (11,"(2i,1d15.7)") i, hist(i), weight(i)
!				write (11,"(2i,1d15.7)") i, hist(i), weight(i)
			enddo
			close(11)
		end subroutine cluster_distribution


		subroutine compute_force_per_cluster
			implicit none
			integer :: i, j, idim

			do i=1, ncluster
				cluster(i)%force(:) = zero
				do j=1, cluster(i)%num
					do idim=1, ndim
						cluster(i)%force(idim) = cluster(i)%force(idim) + force(cluster(i)%list(j),idim)
					enddo
				enddo
				do idim=1, ndim
					cluster(i)%force(idim) = cluster(i)%force(idim) / cluster(i)%num
				enddo
			enddo
		end subroutine compute_force_per_cluster

	end subroutine find_clusters
end module clustermod


#if 0
		subroutine projected_area
			implicit none
			real(prcn) :: dist(2), sep, alpha
			integer :: icluster, idim, jdim, i, j, dim(2)
			real(prcn), allocatable :: cluster2d(:,:)

			do icluster=1, ncluster
				cluster(icluster)%area(:) = zero
!				cord(idim) = 0

				allocate(cluster2d(cluster(icluster)%num,2))
				cluster2d = zero

				do idim=1, ndim
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
							cluster2d(i,jdim) = xc( cluster(icluster)%list(i),dim(jdim) ) / dbydx
						enddo
					enddo

					cluster(icluster)%area(idim) = pi*dia_phys**2 * cluster(icluster)%num

					do i=1, cluster(icluster)%num-1
						do j=i+1, cluster(icluster)%num
							do jdim=1, 2
								dist(jdim) = abs( cluster2d(j,jdim)-cluster2d(i,jdim) )

!								if (jdim==1 .and. dist(jdim)>cord(idim)) cord(idim) = dist(jdim)

								if (dist(jdim)>rmax) dist(jdim) = 2*rmax-dist(jdim)
							enddo

							sep = sqrt( dot_product(dist,dist) )
							
							if (sep<dia_phys) then
								alpha = asin(sep / dia_phys)

								cluster(icluster)%area(idim) = cluster(icluster)%area(idim) - 0.5 * (dia_phys/2)**2 * (2*alpha-sin(2*alpha))
							endif
						enddo
					enddo
				enddo
				deallocate(cluster2d)

!				write (600,"(12d15.7)") surface(:)/(pi*rad_tmp(1)**2), surface(1)/surface(2), surface(1)/surface(3), surface(2)/surface(3), cord(:)/(2*rad_tmp(1)), cord(1)/cord(2), cord(1)/cord(3), cord(2)/cord(3)
			enddo
		end subroutine projected_area
#endif
