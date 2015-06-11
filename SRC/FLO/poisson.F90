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


module poisson
#include "ibm.h"
	USE precision 
	USE constants 
	USE global_data
	USE bcsetarrays
	use fftw3_interface

	IMPLICIT NONE 
	!///////////////////////////////////////////////////////////////////
	!	Calculate the pressure field using AB. The boundary conditions
	!	are such that only the flow field from the previous two
	!	timesteps (at the boundaries) is needed to solve the system.
	!	The pressure is estimated at (t+dt)
	!--------------------------------------------------------------------
contains 

	subroutine pressure
		implicit none 
		!-----------------------------------------------------------------------
		!	local variables
		complex(prcn) :: tmpc
		real(prcn) :: b
		complex(prcn) :: wtmp
		integer :: i, j, k, idim

#if !PARALLEL
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					if (w2(i,j,k)>small_number) then
						b = -w2(i,j,k)
						tmpc = czero
						do idim=1, ndim
							if (idim == 1) then
								wtmp = wx(i)
							elseif (idim == 2) then
								wtmp = wy(j)
							elseif (idim == 3) then
								wtmp = wz(k)
							endif
#else
		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					if (w2(i,j,k)>small_number) then
						b = -w2(i,j,k)
						tmpc = czero
						do idim=1, ndim
							if (idim == 1) then
								wtmp = wy(j)
							elseif (idim == 2) then
								wtmp = wz(k)
							elseif (idim == 3) then
								wtmp = wx(i)
							endif
#endif

							tmpc = tmpc + wtmp * (coef(1,3)*nl(i,j,k,idim) + coef(1,4)*onl(i,j,k,idim) + ff(i,j,k,idim))
						enddo
						tmpc = tmpc / b
						p(i,j,k) = p(i,j,k) + prf * (tmpc-p(i,j,k))
					endif
				enddo
			enddo
		enddo
	end subroutine pressure

	subroutine divcorr
		implicit none 
		!-----------------------------------------------------------------------
		!	local variables
		real(prcn) ::  b
		complex(prcn) :: tmpc, wtmp, div
		integer :: i,j,k, idim

		!new gauss_phi= .true.
#if !PARALLEL
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					if (w2(i,j,k)>small_number) then
						b = -dt*w2(i,j,k)
						tmpc = czero
						do idim=1, ndim
							if (idim == 1) then
								wtmp = wx(i)
							elseif (idim == 2) then
								wtmp = wy(j)
							elseif (idim == 3) then
								wtmp = wz(k)
							endif
#else
		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					if (w2(i,j,k)>small_number) then
						b = -dt*w2(i,j,k)
						tmpc = czero
						do idim=1, ndim
							if (idim == 1) then
								wtmp = wy(j)
							elseif (idim == 2) then
								wtmp = wz(k)
							elseif (idim == 3) then
								wtmp = wx(i)
							endif
#endif
							!-----------------------------------------------------
							!	calculate divergence of intermediate velocity field
							!	form divergence in Fourier space on pressure grid
							tmpc = tmpc + wtmp * u(i,j,k,idim)
						enddo

						!---------------------------------
						!    velocity corection
						div = tmpc/b
						do idim=1, ndim
#if !PARALLEL
							if (idim == 1) then
								wtmp = wx(i)
							elseif (idim == 2) then
								wtmp = wy(j)
							elseif (idim == 3) then
								wtmp = wz(k)
							endif
#else
							if (idim == 1) then
								wtmp = wy(j)
							elseif (idim == 2) then
								wtmp = wz(k)
							elseif (idim == 3) then
								wtmp = wx(i)
							endif
#endif
							u(i,j,k,idim) = u(i,j,k,idim) - dt*wtmp*div
						enddo

						!---------------------------------
						!    pressure corection
						p(i,j,k) = p(i,j,k) + div + half*dt*vis*w2(i,j,k)*div
					endif
				enddo
			enddo
		enddo
		!new gauss_phi= .false.
	end subroutine divcorr
end module poisson
