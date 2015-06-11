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


module nl_allflow
#include "ibm.h" 
	use precision 
	use constants 
	use global_data
	use machine
	use field_tmp_arrays
	use fftw3_interface
	use nlmainarrays, Only : ubcp
	implicit none 
	private
	public :: form_nl
	real(prcn) ::  u_max, v_max, w_max
	real(prcn) ::  umax_loc, vmax_loc, wmax_loc
contains 
	subroutine form_nl
		implicit none 
		integer :: i,j,k, dim1, dim2, ii, jj, kk
		real(prcn) :: slope_factor, mesh_veltemp(ndim) !max_speed, Dplus_vj, Dminus_vj, phi_ri, 

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
end module nl_allflow

