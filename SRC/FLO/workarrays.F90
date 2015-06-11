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

Module nlarrays 
	Use precision 
	Use constants 
	Use global_data
	implicit none

#if PARALLEL
	!REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: uatminus1, uatnxp2
#endif
end Module nlarrays

Module bcsetarrays
	Use precision 
	Use constants 
	Use global_data
	implicit none

	Real(prcn),  DIMENSION(:,:,:,:), ALLOCATABLE, target :: omega, fr, ppr, diffn
end Module bcsetarrays

Module nlmainarrays
	Use precision 
	Use constants 
	Use global_data
	implicit none

	real(prcn), DIMENSION(:,:,:,:), ALLOCATABLE,  Target ::  ubc, nlbc, onlbc
	real(prcn), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pbc

	real(prcn), Dimension(:,:,:,:), pointer ::  onlbcp, nlbcp, ubcp
	real(prcn), Dimension(:,:,:), pointer ::  pbcp
end Module nlmainarrays

module field_tmp_arrays
	use precision
	implicit none

	real(prcn), allocatable :: urtmp(:,:,:)
	complex(prcn), allocatable :: uftmp(:,:,:) !, uftmp2(:,:,:), uftmp3(:,:,:)
end module field_tmp_arrays
	



