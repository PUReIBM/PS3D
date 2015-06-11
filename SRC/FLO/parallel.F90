!    PUReIBM-PS3D is a three-dimensional psudeo-spectral particle-resolved
!    direct numerical simulation solver for detailed analysis of homogeneous
!    fixed and freely evolving fluid-particle suspensions. PUReRIBM-PS3D
!    is a continuum Navier-Stokes solver based on Cartesian grid that utilizes
!    Immeresed Boundary method to represent particle surfuces.
!    Copyright (C) 2015, Shankar Subramaniam, Rahul Garg, Sudheer Tenneti, Bo Sun, Mohammad Mehrabadi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    For acknowledgement, please refer to the following publications:
!    (1) TENNETI, S. & SUBRAMANIAM, S., 2014, Particle-resolved direct numerical
!        simulation for gas–solid flow model development. Annu. Rev. Fluid Mech.
!        46 (1), 199–230.
!    (2) SUBRAMANIAM, S., MEHRABADI, M., HORWITZ, J. & MANI, A., 2014, Developing
!        improved Lagrangian point particle models of gas–solid flow from
!        particle-resolved direct numerical simulation. In Studying Turbulence
!        Using Numerical Simulation Databases-XV, Proceedings of the CTR 2014
!        Summer Program, pp. 5–14. Center for Turbulence Research, Stanford
!        University, CA.


module parallel
#include "ibm.h"
	use global_data
	use fftw3_interface
	implicit none

contains

	subroutine initialize_parallel
		implicit none
		integer :: myid_old

#if PARALLEL
		communicator_done = .false.
		PARALLEL_START()
		call pfft_init()

		GET_NPROCS(comm_group, nproc)
		GET_PROCESSOR_RANK(comm_group, myid)

		allocate(nproc_array(2))
		nproc_array(1) = nprocz
		nproc_array(2) = nprocy
		
		CREATE_CART_TOPOLOGY(comm_group, ndim-1, nproc_array, xperiodic, .true., decomp_group)  
		GET_PROCESSOR_RANK(decomp_group, myid_old)

		err_code = pfft_create_procmesh_2d(decomp_group, nprocz, nprocy, comm_cart_2d)
		GET_PROCESSOR_RANK(comm_cart_2d, myid)

		if (myid/=myid_old) then
			call screen_separator(60,'!')
			write (*,"(2(1a,1i))") "WARNING, MYID_OLD = ", myid_old, ", MYID = ", myid
			call screen_separator(60,'-')
		endif

		if (err_code/=0 .or. nprocy*nprocz/=nproc) then
			if (I_AM_NODE_ZERO) then
				write(*,"(1a)") "ERROR: THIS JOB FAILED TO INITIALIZE BECAUSE:"
				write(*,"(1a,1i)") "  NPROC  = ", nproc
				write(*,"(1a,1i)") "  NPROCX = ", nprocy
				write(*,"(1a,1i)") "  NPROCY = ", nprocz
			endif
			PARALLEL_FINISH()
			stop
		endif

		if (I_AM_NODE_ZERO) call screen_separator(60,'P')
		BARRIER(comm_cart_2d)

		GET_SHIFT_PROCS(comm_cart_2d,1,1,fromy,toy)
		GET_SHIFT_PROCS(comm_cart_2d,0,1,fromz,toz)
		
		if (I_AM_NODE_ZERO) write (*,*) "RANKS AND NEIGHBORS..."
		write (*,"(5(1a10,1i6))") "MYID = ",  myid, ", FORMY = ", fromy, ", TOY = ", toy, ", FROMZ = ", fromz, ", TOZ = ", toz
		BARRIER(comm_cart_2d)

		if (I_AM_NODE_ZERO) call screen_separator(60,'-')
		communicator_done = .true.
#else
		nproc = 1
		nprocy = 1
		nprocz = 1
		myid = 0
#endif
	end subroutine initialize_parallel


	function point_in_this_domain(index_in,index_out,position_in,position_out)
		implicit none
		logical :: point_in_this_domain
		integer , intent(in)  :: index_in(ndim)
		integer , intent(out) :: index_out(ndim)
		real(prcn) , intent(in)  :: position_in(ndim)
		real(prcn) , intent(out) :: position_out(ndim)

		integer :: idim

		point_in_this_domain = .true.

		do idim=1, ndim
			if (index_in(idim) < 1) then 
				index_out(idim)    = index_in(idim)    + global_n(idim)
				position_out(idim) = position_in(idim) + global_n(idim)
			elseif (index_in(idim) > global_n(idim)) then 
				index_out(idim)    = index_in(idim)    - global_n(idim)
				position_out(idim) = position_in(idim) - global_n(idim)
			else
				index_out(idim)    = index_in(idim)
				position_out(idim) = position_in(idim)
			endif
		enddo

		do idim=1, ndim
			if (local_i_start(idim)<=index_out(idim) .and. index_out(idim)<=local_i_end(idim)) then

			else
				point_in_this_domain = .false.
				exit
			endif
		enddo

		if (point_in_this_domain) then
			index_out(:)    = index_out(:)    - local_i_start(:) + 1
			position_out(:) = position_out(:) - local_i_start(:) + 1
		endif
	end function point_in_this_domain


	function point_in_this_ghost_domain(index_in,index_out,position_in,position_out)
		implicit none
		logical :: point_in_this_ghost_domain
		integer , intent(in)  :: index_in(ndim)
		integer , intent(out) :: index_out(ndim)
		real(prcn) , intent(in)  :: position_in(ndim)
		real(prcn) , intent(out) :: position_out(ndim)

		integer :: index_tmp(ndim)
		real(prcn) :: position_tmp(ndim)
		integer :: idim

		point_in_this_ghost_domain = .true.

		do idim=1, ndim
			if (index_in(idim) <= 1) then 
				index_tmp(idim)    = index_in(idim)    + global_n(idim)
				position_tmp(idim) = position_in(idim) + global_n(idim)
			elseif (index_in(idim) >= global_n(idim)-1) then 
				index_tmp(idim)    = index_in(idim)    - global_n(idim)
				position_tmp(idim) = position_in(idim) - global_n(idim)
			else
				index_tmp(idim)    = index_in(idim)
				position_tmp(idim) = position_in(idim)
			endif
		enddo

		do idim=1, ndim
			if (local_i_start(idim)-2<=index_in(idim) .and. index_in(idim)<=local_i_end(idim)+1) then
				index_out(idim)    = index_in(idim)    - local_i_start(idim) + 1
				position_out(idim) = position_in(idim) - local_i_start(idim) + 1
			elseif (local_i_start(idim)-2<=index_tmp(idim) .and. index_tmp(idim)<=local_i_end(idim)+1) then
				index_out(idim)    = index_tmp(idim)    - local_i_start(idim) + 1
				position_out(idim) = position_tmp(idim) - local_i_start(idim) + 1
			else
				point_in_this_ghost_domain = .false.
				exit
			endif
		enddo

#if !PARALLEL
		if (.not.point_in_this_ghost_domain) then
			write (*,"(1A,3i8,3d15.7)") "A NODE NOT IN THIS DOMAIN : ", index_in, position_in
		endif
#endif
		! THIS POINT SHOULD ONLY BE IN A SUB-DOMAIN.
		! CAN BE CHECKED BY A GLOBAL SUM. THE RESULT SHOULD BE 1
	end function point_in_this_ghost_domain


	function cell_in_this_domain(index_in,index_out,position_in,position_out)
		implicit none
		logical :: cell_in_this_domain
		integer , intent(in)  :: index_in(ndim)
		integer , intent(out) :: index_out(ndim)
		real(prcn) , intent(in)  :: position_in(ndim)
		real(prcn) , intent(out) :: position_out(ndim)

		integer :: index_tmp(ndim)
		real(prcn) :: position_tmp(ndim)
		integer :: idim

		cell_in_this_domain = .true.

		do idim=1, ndim
			if (index_in(idim) < 1) then 
				index_tmp(idim)    = index_in(idim)    + global_n(idim)
				position_tmp(idim) = position_in(idim) + global_n(idim)
			elseif (index_in(idim) >= global_n(idim)) then 
				index_tmp(idim)    = index_in(idim)    - global_n(idim)
				position_tmp(idim) = position_in(idim) - global_n(idim)
			else
				index_tmp(idim)    = index_in(idim)
				position_tmp(idim) = position_in(idim)
			endif
		enddo

#if 0
if (check_point) then
	write (*,"(1a,3i)") "index_in  = ", index_in
	write (*,"(1a,3i)") "index_tmp = ", index_tmp

	write (*,"(1a,3d15.7)") "position_in  = ", position_in
	write (*,"(1a,3d15.7)") "position_tmp = ", position_tmp
endif
#endif

		do idim=1, ndim
			if (local_i_start(idim)-1<=index_in(idim) .and. index_in(idim)<=local_i_end(idim)) then
				index_out(idim)    = index_in(idim)    - local_i_start(idim) + 1
				position_out(idim) = position_in(idim) - local_i_start(idim) + 1

#if 0
if (check_point) then
	write (*,"(1a,5i)") "index_out(idim)  = ", 1, idim, index_out(idim)
endif
#endif

			elseif (local_i_start(idim)-1<=index_tmp(idim) .and. index_tmp(idim)<=local_i_end(idim)) then
				index_out(idim)    = index_tmp(idim)    - local_i_start(idim) + 1
				position_out(idim) = position_tmp(idim) - local_i_start(idim) + 1

#if 0
if (check_point) then
	write (*,"(1a,5i)") "index_out(idim)  = ", 2, idim, index_out(idim)
endif
#endif

			else
				cell_in_this_domain = .false.

				index_out = -10
				position_out = -10

				exit
			endif
		enddo

#if 0
if (check_point) then
	write (*,*) "cell_in_this_domain = ", cell_in_this_domain
	write (*,"(1a,3i)") "index_in  = ", index_in
	write (*,"(1a,3i)") "index_out = ", index_out

	write (*,"(1a,3d15.7)") "position_in  = ", position_in
	write (*,"(1a,3d15.7)") "position_tmp = ", position_out
endif
#endif


	end function cell_in_this_domain






#if 0
	function cell_in_this_ghost_domain(index_in,index_out,position_in,position_out)
		implicit none
		logical :: cell_in_this_ghost_domain
		integer , intent(in)  :: index_in(ndim)
		integer , intent(out) :: index_out(ndim)
		real(prcn) , intent(in)  :: position_in(ndim)
		real(prcn) , intent(out) :: position_out(ndim)

		integer :: index_tmp(ndim)
		real(prcn) :: position_tmp(ndim)
		integer :: idim

		cell_in_this_ghost_domain = .true.

		do idim=1, ndim
			if (index_in(idim) <= 1) then 
				index_tmp(idim)    = index_in(idim)    + global_n(idim)
				position_tmp(idim) = position_in(idim) + global_n(idim)
			elseif (index_in(idim) >= global_n(idim)-1) then 
				index_tmp(idim)    = index_in(idim)    - global_n(idim)
				position_tmp(idim) = position_in(idim) - global_n(idim)
			else
				index_tmp(idim)    = index_in(idim)
				position_tmp(idim) = position_in(idim)
			endif
		enddo

		do idim=1, ndim

			if (local_i_start(idim)-2<=index_tmp(idim) .and. index_tmp(idim)<=local_i_end(idim)+1) then
				index_out(idim)    = index_in(idim)    - local_i_start(idim) + 1
				position_out(idim) = position_in(idim) - local_i_start(idim) + 1
			elseif (local_i_start(idim)-2<=index_in(idim) .and. index_in(idim)<=local_i_end(idim)+1) then
				index_out(idim)    = index_tmp(idim)    - local_i_start(idim) + 1
				position_out(idim) = position_tmp(idim) - local_i_start(idim) + 1
			else
				cell_in_this_ghost_domain = .false.
				exit
			endif
		enddo

#if !PARALLEL
		if (.not.cell_in_this_ghost_domain) then
			write (*,"(1A,3i8,3d15.7)") "A NODE NOT IN THIS DOMAIN : ", index_in, position_in
		endif
#endif
		! THIS POINT SHOULD ONLY BE IN A SUB-DOMAIN.
		! CAN BE CHECKED BY A GLOBAL SUM. THE RESULT SHOULD BE 1
	end function cell_in_this_ghost_domain
#endif


	function in_how_many_domains(log_in)
		implicit none
		logical, intent(in) :: log_in
		integer :: in_how_many_domains
		integer :: count

		if (log_in) then
			count = 1
		else
			count = 0
		endif

		GLOBAL_INT_SUM(count, in_how_many_domains, 1, comm_cart_2d)
	end function in_how_many_domains


	subroutine communicate_fluid_atijk
		implicit none

		if (nprocz>1) then
			!COMMUNICATION ALONG -Z
			VECSENDRECV(fluid_atijk(0,0,1), 1, xy_logical, fromz, 1, &
						&	fluid_atijk(0,0,local_ni(3)+1), 1, toz, 1, comm_cart_2d, status)
			!COMMUNICATION ALONG +Z
			VECSENDRECV(fluid_atijk(0,0,local_ni(3)), 1, xy_logical, toz, 2, &
						&	fluid_atijk(0,0,0), 1, fromz, 2, comm_cart_2d, status)
		else
			!COMMUNICATION ALONG Z
			fluid_atijk(:,:,local_ni(3)+1) = fluid_atijk(:,:,1)
			fluid_atijk(:,:,0) = fluid_atijk(:,:,local_ni(3))
		endif

		if (nprocy>1) then
			!COMMUNICATION ALONG -Y
			VECSENDRECV(fluid_atijk(0,1,0), 1, xz_logical, fromy, 3, &
						&	fluid_atijk(0,local_ni(2)+1,0), 1, toy, 3, comm_cart_2d, status)
			!COMMUNICATION ALONG +Y
			VECSENDRECV(fluid_atijk(0,local_ni(2),0), 1, xz_logical, toy, 4, &
						&	fluid_atijk(0,0,0), 1, fromy, 4, comm_cart_2d, status)
		else
			!COMMUNICATION ALONG Y
			fluid_atijk(:,local_ni(2)+1,:) = fluid_atijk(:,1,:)
			fluid_atijk(:,0,:) = fluid_atijk(:,local_ni(2),:)
		endif

		!COMMUNICATION ALONG X
		fluid_atijk(local_ni(1)+1,:,:) = fluid_atijk(1,:,:)
		fluid_atijk(0,:,:) = fluid_atijk(local_ni(1),:,:)

	end subroutine communicate_fluid_atijk


	subroutine communicate_velocity(ur_in)
		implicit none
		real(prcn), intent(inout) :: ur_in(-1:local_ni(1)+2,-1:local_ni(2)+2,-1:local_ni(3)+2,1:ndim)
		integer :: idim

		do idim=1, ndim
			if (nprocz>1) then
				!COMMUNICATION ALONG -Z, FIRST LAYER
				VECSENDRECV(ur_in(-1,-1,1,idim), 1, xy_vel, fromz, 5, &
							&	ur_in(-1,-1,local_ni(3)+1,idim), 1, toz, 5, comm_cart_2d, status)
				!COMMUNICATION ALONG +Z, FIRST LAYER
				VECSENDRECV(ur_in(-1,-1,local_ni(3),idim), 1, xy_vel, toz, 6, &
							&	ur_in(-1,-1,0,idim), 1, fromz, 6, comm_cart_2d, status)

				!COMMUNICATION ALONG -Z, SECOND LAYER
				VECSENDRECV(ur_in(-1,-1,2,idim), 1, xy_vel, fromz, 7, &
							&	ur_in(-1,-1,local_ni(3)+2,idim), 1, toz, 7, comm_cart_2d, status)
				!COMMUNICATION ALONG +Z, SECOND LAYER
				VECSENDRECV(ur_in(-1,-1,local_ni(3)-1,idim), 1, xy_vel, toz, 8, &
							&	ur_in(-1,-1,-1,idim), 1, fromz, 8, comm_cart_2d, status)
			else
				!COMMUNICATION ALONG Z, FIRST LAYER
				ur_in(:,:,0,idim) = ur_in(:,:,local_ni(3),idim)
				ur_in(:,:,local_ni(3)+1,idim) = ur_in(:,:,1,idim)

				!COMMUNICATION ALONG Z, SECOND LAYER
				ur_in(:,:,-1,idim) = ur_in(:,:,local_ni(3)-1,idim)
				ur_in(:,:,local_ni(3)+2,idim) = ur_in(:,:,2,idim)
			endif

			if (nprocy>1) then
				!COMMUNICATION ALONG -Y, FIRST LAYER
				VECSENDRECV(ur_in(-1,1,-1,idim), 1, xz_vel, fromy, 9, &
							&	ur_in(-1,local_ni(2)+1,-1,idim), 1, toy, 9, comm_cart_2d, status)
				!COMMUNICATION ALONG +Y, FIRST LAYER
				VECSENDRECV(ur_in(-1,local_ni(2),-1,idim), 1, xz_vel, toy, 10, &
							&	ur_in(-1,0,-1,idim), 1, fromy, 10, comm_cart_2d, status)

				!COMMUNICATION ALONG -Y, SECOND LAYER
				VECSENDRECV(ur_in(-1,2,-1,idim), 1, xz_vel, fromy, 11, &
							&	ur_in(-1,local_ni(2)+2,-1,idim), 1, toy, 11, comm_cart_2d, status)
				!COMMUNICATION ALONG +Y, SECOND LAYER
				VECSENDRECV(ur_in(-1,local_ni(2)-1,-1,idim), 1, xz_vel, toy, 12, &
							&	ur_in(-1,-1,-1,idim), 1, fromy, 12, comm_cart_2d, status)
			else
				!COMMUNICATION ALONG Y, FIRST LAYER
				ur_in(:,0,:,idim) = ur_in(:,local_ni(2),:,idim)
				ur_in(:,local_ni(2)+1,:,idim) = ur_in(:,1,:,idim)

				!COMMUNICATION ALONG Y, SECOND LAYER
				ur_in(:,-1,:,idim) = ur_in(:,local_ni(2)-1,:,idim)
				ur_in(:,local_ni(2)+2,:,idim) = ur_in(:,2,:,idim)
			endif

			!COMMUNICATION ALONG X, FIRST LAYER
			ur_in(0,:,:,idim) = ur_in(local_ni(1),:,:,idim)
			ur_in(local_ni(1)+1,:,:,idim) = ur_in(1,:,:,idim)

			!COMMUNICATION ALONG X, SECOND LAYER
			ur_in(-1,:,:,idim) = ur_in(local_ni(1)-1,:,:,idim)
			ur_in(local_ni(1)+2,:,:,idim) = ur_in(2,:,:,idim)
		enddo
	end subroutine communicate_velocity


	subroutine communicate_in_gohst_domain(input)
		implicit none
		real(prcn), intent(inout) :: input(0:local_ni(1)+1,0:local_ni(2)+1,0:local_ni(3)+1)

		if (nprocz>1) then
			!COMMUNICATION ALONG -Z
			VECSENDRECV(input(0,0,1), 1, xy_double, fromz, 13, &
						&	input(0,0,local_ni(3)+1), 1, toz, 13, comm_cart_2d, status)
			!COMMUNICATION ALONG +Z
			VECSENDRECV(input(0,0,local_ni(3)), 1, xy_double, toz, 14, &
						&	input(0,0,0), 1, fromz, 14, comm_cart_2d, status)
		else
			!COMMUNICATION ALONG Z
			input(:,:,0) = input(:,:,local_ni(3))
			input(:,:,local_ni(3)+1) = input(:,:,1)
		endif

		if (nprocy>1) then
			!COMMUNICATION ALONG -Y
			VECSENDRECV(input(0,1,0), 1, xz_double, fromy, 15, &
						&	input(0,local_ni(2)+1,0), 1, toy, 15, comm_cart_2d, status)
			!COMMUNICATION ALONG +Y
			VECSENDRECV(input(0,local_ni(2),0), 1, xz_double, toy, 16, &
						&	input(0,0,0), 1, fromy, 16, comm_cart_2d, status)
		else
			!COMMUNICATION ALONG Y
			input(:,0,:) = input(:,local_ni(2),:)
			input(:,local_ni(2)+1,:) = input(:,1,:)
		endif

		!COMMUNICATION ALONG X
		input(0,:,:) = input(local_ni(1),:,:)
		input(local_ni(1)+1,:,:) = input(1,:,:)
	end subroutine communicate_in_gohst_domain


	subroutine communicate_ib_forces
		use bcsetarrays, only : fr, diffn
		implicit none
		integer :: idim

		do idim=1, ndim
			if (nprocz>1) then
				!COMMUNICATION ALONG -Z
				VECSENDRECV(fr(0,0,0,idim), 1, xy_double, fromz, 17, &
							&	diffn(0,0,local_ni(3)+1,idim), 1, toz, 17, comm_cart_2d, status)
				!COMMUNICATION ALONG +Z
				VECSENDRECV(fr(0,0,local_ni(3)+1,idim), 1, xy_double, toz, 18, &
							&	diffn(0,0,0,idim), 1, fromz, 18, comm_cart_2d, status)

				fr(:,:,local_ni(3),idim) = fr(:,:,local_ni(3),idim) + diffn(:,:,local_ni(3)+1,idim)
				fr(:,:,1,idim) = fr(:,:,1,idim) + diffn(:,:,0,idim)


			else
				!COMMUNICATION ALONG Z
				fr(:,:,local_ni(3),idim) = fr(:,:,local_ni(3),idim) + fr(:,:,0,idim)
				fr(:,:,1,idim)           = fr(:,:,1,idim)           + fr(:,:,local_ni(3)+1,idim)
			endif

			if (nprocy>1) then
				!COMMUNICATION ALONG -Y
				VECSENDRECV(fr(0,0,0,idim), 1, xz_double, fromy, 19, &
							&	diffn(0,local_ni(2)+1,0,idim), 1, toy, 19, comm_cart_2d, status)
				!COMMUNICATION ALONG +Y
				VECSENDRECV(fr(0,local_ni(2)+1,0,idim), 1, xz_double, toy, 20, &
							&	diffn(0,0,0,idim), 1, fromy, 20, comm_cart_2d, status)

				fr(:,local_ni(2),:,idim) = fr(:,local_ni(2),:,idim) + diffn(:,local_ni(2)+1,:,idim)
				fr(:,1,:,idim) = fr(:,1,:,idim) + diffn(:,0,:,idim)
			else
				!COMMUNICATION ALONG Y
				fr(:,local_ni(2),:,idim) = fr(:,local_ni(2),:,idim) + fr(:,0,:,idim)
				fr(:,1,:,idim)           = fr(:,1,:,idim)           + fr(:,local_ni(2)+1,:,idim)
			endif

			!COMMUNICATION ALONG X
			fr(local_ni(1),:,:,idim) = fr(local_ni(1),:,:,idim) + fr(0,:,:,idim)
			fr(1,:,:,idim)           = fr(1,:,:,idim)           + fr(local_ni(1)+1,:,:,idim)
		enddo

	end subroutine communicate_ib_forces
end module parallel

