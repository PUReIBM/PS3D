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

MODULE usteptime
#include "ibm.h"
	USE precision 
	USE constants 
	USE global_data
	use fftw3_interface
	implicit none 
Contains
	SUBROUTINE ustep
		IMPLICIT NONE 
		!-----------------------------------------------------------------------
		!	local variables
		REAL(prcn) :: b
		COMPLEX(prcn) ::  r, wtmp
		INTEGER :: i, j, k, idim

		!gauss_u = .true.
		do idim=1, ndim
#if !PARALLEL
			do k=1, local_no(3)
				do j=1, local_no(2)
					do i=1, local_no(1)
						if (w2(i,j,k)>small_number) then
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
						if(w2(i,j,k)>small_number) then
							if (idim == 1) then
								wtmp = wy(j)
							elseif (idim == 2) then
								wtmp = wz(k)
							elseif (idim == 3) then
								wtmp = wx(i)
							endif
#endif
							!-----------------------------------------------------------------------
							!	the left-hand-side
							b = one + half*dt*vis*w2(i,j,k)

							!-----------------------------------------------------------------------
							!	diffusion
							r = u(i,j,k,idim) * (one - half*dt*vis*w2(i,j,k))

							!-----------------------------------------------------------------------
							!	include convection terms
							!       sign change infront of dt done..
							!       to be consistent with Jamals convention
							!       remember nl term is negative..
							!       so negative sign, which should have been in front of dt,
							!       is now absorbed into the nl term
							r = r + dt * (coef(1,3)*nl(i,j,k,idim)+coef(1,4)*onl(i,j,k,idim)) + dt*ff(i,j,k,idim)*force_factor

							!-----------------------------------------------------------------------
							!	pressure terms
							r = r - dt * p(i,j,k)*wtmp

							u(i,j,k,idim) = r/b
						endif
					enddo
				enddo
			enddo
		enddo
		!new gauss_u = .false.
	END SUBROUTINE ustep
END MODULE usteptime
