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


module fftw3_interface
#include "ibm.h"
	use, intrinsic :: iso_c_binding
	USE precision
	USE constants
	USE global_data
	use general_funcs, only : screen_separator
	IMPLICIT NONE
#if !PARALLEL
	include "/home/mehr/FFTW3/include/fftw3.f03"
#else
	include "/home/mehr/FFTW3/include/fftw3-mpi.f03"	
	include "/home/mehr/PFFT/include/pfft.f03"
	integer(C_SIZE_T) :: local_ni_fftw(ndim)
#endif

	type(C_PTR) :: real_ptr, complex_ptr, plan_r2c, plan_c2r
	integer(C_INTPTR_T) :: real_size, complex_size
	real(C_DOUBLE),            pointer :: real_array(:,:,:)
	complex(C_DOUBLE_COMPLEX), pointer :: complex_array(:,:,:)

	integer(C_SIZE_T) :: global_n(ndim), global_n_c(ndim)
	integer(C_SIZE_T) :: local_ni(ndim), local_i_start(ndim), local_i_end(ndim), local_ni_c(ndim), local_i_start_c(ndim)
	integer(C_SIZE_T) :: local_no(ndim), local_o_start(ndim), local_o_end(ndim), local_no_c(ndim), local_o_start_c(ndim)

	integer, allocatable :: xlocals(:), ylocals(:), zlocals(:)
	integer, allocatable :: xstarts(:), ystarts(:), zstarts(:)
	integer, allocatable :: xends(:), yends(:), zends(:)

	integer, allocatable :: ixlocals(:), iylocals(:), izlocals(:)
	integer, allocatable :: ixstarts(:), iystarts(:), izstarts(:)
	integer, allocatable :: ixends(:), iyends(:), izends(:)
	integer :: ierr

contains

	subroutine initiate_fftw3
		use field_tmp_arrays
		implicit none
		integer :: i, j

		global_n(1) = mx
		global_n(2) = my
		global_n(3) = mz

		global_n_c(1) = global_n(3)
		global_n_c(2) = global_n(2)
		global_n_c(3) = global_n(1)

#if !PARALLEL
		!LENGTH OF REAL SPACE DIMENSIONS
		local_ni(:) = global_n(:)
		local_i_start(:) = 1
		local_i_end(:) = local_i_start(:)+local_ni(:)-1

		local_ni_c(1) = global_n_c(1)
		local_ni_c(2) = global_n_c(2)
		local_ni_c(3) = global_n_c(3)

		!LENGTH OF FOURIER SPACE DIMENSIONS
		local_no(1) = global_n(1)/2+1
		local_no(2) = global_n(2)
		local_no(3) = global_n(3)

		local_o_start (:) = 1
		local_o_end(:) = local_o_start(:)+local_no(:)-1

		complex_size = local_no(1)*local_no(2)*local_no(3)
		real_size    = local_ni(1)*local_ni(2)*local_ni(3)

		real_ptr    = fftw_alloc_real(real_size)
		complex_ptr = fftw_alloc_complex(complex_size)


		call c_f_pointer(real_ptr,       real_array, [local_ni(1),local_ni(2),local_ni(3)])
		call c_f_pointer(complex_ptr, complex_array, [local_no(1),local_no(2),local_no(3)])

		plan_r2c = fftw_plan_dft_r2c_3d(local_ni_c(1),local_ni_c(2),local_ni_c(3), real_array, complex_array, FFTW_MEASURE+FFTW_DESTROY_INPUT)
		plan_c2r = fftw_plan_dft_c2r_3d(local_ni_c(1),local_ni_c(2),local_ni_c(3), complex_array, real_array, FFTW_MEASURE+FFTW_DESTROY_INPUT)

#else
		complex_size = pfft_local_size_dft_r2c_3d(global_n_c, comm_cart_2d, PFFT_TRANSPOSED_OUT, &
							& local_ni_c, local_i_start_c, local_no_c, local_o_start_c)
		real_size    = 2*complex_size

		!NOTE THAT IN PARALLEL MODE THE OUTPUT OF TRANSFORMATION TO FOURIER SPACE IS TRANSPOSED.
		!ANY 3D ARRAY IS TRACABLE AS FOLLOWS:
		!NX x NY/PY x NZ/PZ ---> (FORTRAN TO C) NZ/PZ x NY/PY x NX ---> (PFFT REFER TO THE MANUAL)
		!MY/PZ x MX/PY x MZ ---> (C TO FORTRAN) MZ x MX/PY x MY/PZ
		!THIS MEAN THAT IN ARRAYS IN FOURIER SPACE SHOULD BE ALLOCATED AND TREATED ACCORDING TO
		!THE DATA STRUCTURE MENTIONED ABOVE

		!IN FOURIER TRANSFORM AND EXTRA PADDING IN THE REAL SPACE IS NEEDED FOR INTERNAL OPERATIONS.
		!local_ni_fftw IS DEFINED TO ACCOUNT FOR THIS EXTRA PADING LOCAL_NI(1) = 2*(GLOBAL_NI(1)/2+1).
		!HOWEVER, FLOW ARRAYS SHOULD BE ALLOCATED EXACTLY AS THE REAL DOMAIN+GOHST CELLS.
		!THEREFOR LOCAL_NI IS ALSO DEFINED FOR THIS PURPOSE, WHICH IS LOCAL_NI(1) = GLOBALE_N(1)

		local_ni_fftw(1) = local_ni_c(3)
		local_ni_fftw(2) = local_ni_c(2)
		local_ni_fftw(3) = local_ni_c(1)

		local_ni(:) = local_ni_fftw(:)
		if (local_ni(1)>global_n(1)) local_ni(1)=global_n(1)

		local_i_start(1) = local_i_start_c(3)+1
		local_i_start(2) = local_i_start_c(2)+1
		local_i_start(3) = local_i_start_c(1)+1

		local_i_end(:) = local_i_start(:)+local_ni(:)-1

		local_no(1) = local_no_c(3)
		local_no(2) = local_no_c(2)
		local_no(3) = local_no_c(1)

		local_o_start(1) = local_o_start_c(3)+1
		local_o_start(2) = local_o_start_c(2)+1
		local_o_start(3) = local_o_start_c(1)+1

		local_o_end(:) = local_o_start(:)+local_no(:)-1

		real_ptr    = pfft_alloc_real(real_size)
		complex_ptr = pfft_alloc_complex(complex_size)

		call c_f_pointer(real_ptr,       real_array, [local_ni_fftw(1),local_ni_fftw(2),local_ni_fftw(3)])
		call c_f_pointer(complex_ptr, complex_array, [local_no(3),local_no(1),local_no(2)])


		plan_r2c = pfft_plan_dft_r2c_3d(global_n_c, real_array, complex_array, comm_cart_2d, PFFT_FORWARD, &
												& PFFT_TRANSPOSED_OUT + PFFT_MEASURE + PFFT_DESTROY_INPUT)

		plan_c2r = pfft_plan_dft_c2r_3d(global_n_c, complex_array, real_array, comm_cart_2d, PFFT_BACKWARD, &
												& PFFT_TRANSPOSED_IN + PFFT_MEASURE + PFFT_DESTROY_INPUT)
#endif

		if (I_AM_NODE_ZERO) then
			call screen_separator(80,'^')
			write (*,*) "LOCAL INDICES"
		endif

		do i=0, nproc-1
			if (myid==i) then
				do j=1, ndim
					write (*,"(3(1i,1a4,2i8,a4))") &
myid,"  | ", local_i_start(j), local_i_end(j), " | ", local_ni(j), " || ", local_o_start(j), local_o_end(j), " | ", local_no(j), " || "

				enddo
				write (*,"(2(1a,1i10))") "REAL_SIZE = ", real_size, ", COMPLEX_SIZE = ", complex_size
				call screen_separator(60,'-')
			endif

			BARRIER(comm_cart_2d)
		enddo

		if (local_no(1)==0.or.local_no(2)==0.or.local_no(3)==0) then
			write (*,*) "WARINIG FROM PROC ", myid,". INEFFICIENT DOMAIN DICOMPISITION IN FFTW2. BETTER CHANGE THE COMPOSITION"
		endif



#if PARALLEL
		allocate(xstarts(0:nproc-1))
		allocate(ystarts(0:nproc-1))
		allocate(zstarts(0:nproc-1))

		allocate(xends(0:nproc-1))
		allocate(yends(0:nproc-1))
		allocate(zends(0:nproc-1))

		call mpi_allgather(local_o_start(1), 1, mpi_int, xstarts, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_o_start(2), 1, mpi_int, ystarts, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_o_start(3), 1, mpi_int, zstarts, 1, mpi_int, comm_cart_2d, ierr)

		call mpi_allgather(local_o_end(1), 1, mpi_int, xends, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_o_end(2), 1, mpi_int, yends, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_o_end(3), 1, mpi_int, zends, 1, mpi_int, comm_cart_2d, ierr)


		allocate(xlocals(0:nproc-1))
		allocate(ylocals(0:nproc-1))
		allocate(zlocals(0:nproc-1))

		call mpi_allgather(local_no(1), 1, mpi_int, xlocals, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_no(2), 1, mpi_int, ylocals, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_no(3), 1, mpi_int, zlocals, 1, mpi_int, comm_cart_2d, ierr)

		!-----------------------------------------------------

		allocate(ixstarts(0:nproc-1))
		allocate(iystarts(0:nproc-1))
		allocate(izstarts(0:nproc-1))

		allocate(ixends(0:nproc-1))
		allocate(iyends(0:nproc-1))
		allocate(izends(0:nproc-1))

		call mpi_allgather(local_i_start(1), 1, mpi_int, ixstarts, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_i_start(2), 1, mpi_int, iystarts, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_i_start(3), 1, mpi_int, izstarts, 1, mpi_int, comm_cart_2d, ierr)

		call mpi_allgather(local_i_end(1), 1, mpi_int, ixends, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_i_end(2), 1, mpi_int, iyends, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_i_end(3), 1, mpi_int, izends, 1, mpi_int, comm_cart_2d, ierr)


		allocate(ixlocals(0:nproc-1))
		allocate(iylocals(0:nproc-1))
		allocate(izlocals(0:nproc-1))

		call mpi_allgather(local_ni(1), 1, mpi_int, ixlocals, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_ni(2), 1, mpi_int, iylocals, 1, mpi_int, comm_cart_2d, ierr)
		call mpi_allgather(local_ni(3), 1, mpi_int, izlocals, 1, mpi_int, comm_cart_2d, ierr)
#endif














		if (I_AM_NODE_ZERO) call screen_separator(80,'H')
	end subroutine initiate_fftw3

	subroutine finilize_fftw3
		implicit none

#if !PARALLEL
		call fftw_destroy_plan(plan_r2c)
		call fftw_destroy_plan(plan_c2r)

		call fftw_free(real_ptr)
		call fftw_free(complex_ptr)
#else
		call pfft_destroy_plan(plan_r2c)
		call pfft_destroy_plan(plan_c2r)
		call pfft_free(real_ptr)
		call pfft_free(complex_ptr)
#endif
	end subroutine finilize_fftw3

	subroutine fftwr2c(real_in, comp_out)
		implicit none
		real(prcn), intent(in), dimension(:,:,:) :: real_in
		complex(prcn), intent(out), dimension(:,:,:) :: comp_out
		integer :: i, j, k

		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					real_array(i,j,k) = real_in(i,j,k)
				enddo
			enddo
		enddo

#if !PARALLEL
		call fftw_execute_dft_r2c(plan_r2c, real_array, complex_array)

		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					comp_out(i,j,k) = complex_array(i,j,k) / (global_n(1)*global_n(2)*global_n(3))
				enddo
			enddo
		enddo
#else
		call pfft_execute(plan_r2c)

		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					comp_out(i,j,k) = complex_array(i,j,k) / (global_n(1)*global_n(2)*global_n(3))
				enddo
			enddo
		enddo
#endif
	end subroutine fftwr2c

	subroutine fftwc2r(comp_in, real_out)
		implicit none
		real(prcn), intent(out), dimension(:,:,:) :: real_out
		complex(prcn), intent(in), dimension(:,:,:) :: comp_in
		integer :: i, j, k

#if !PARALLEL
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					complex_array(i,j,k) = comp_in(i,j,k)
				enddo
			enddo
		enddo

		call fftw_execute_dft_c2r(plan_c2r, complex_array, real_array)
#else
		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					complex_array(i,j,k) = comp_in(i,j,k)
				enddo
			enddo
		enddo

		call pfft_execute(plan_c2r)
#endif

		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					real_out(i,j,k) = real_array(i,j,k)
				enddo
			enddo
		enddo
	end subroutine fftwc2r
end module fftw3_interface
