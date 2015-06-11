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

module steptheflow
#include "ibm.h"
	!///////////////////////////////////////////////////////////////////////
	!	calculate velocity field at next time level
	!	pressure estimation and divergence correction done in place
	!	(for each wavenumber) to save storage
	!-----------------------------------------------------------------------
	USE precision 
	USE constants 
	Use poisson, Only : pressure, divcorr
	use usteptime 
	Use nlmainarrays, Only : pbcp, ubcp
	Use dependent_functions
	USE dem_mod

contains
	subroutine velstep(sflag,rks)
		use init_turb, only : check_divergence, calc_velreal
		use fftw3_interface
		use parallel

		implicit none

		integer :: j, k, i
		integer, intent(in) :: sflag, rks
		complex(prcn) :: usumloc, usum

		CALL pressure
		CALL ustep
		CALL divcorr
		if (debug_check) CALL check_divergence

		! COMMENTING THIS, BECAUSE GRAVITY IS CHANGED IN MPG. IF BOTH GRAVITY AND MPG ARE AVAILABLE, THEN UNCOMMENT THIS
		umean(:) = umean(:) + (-mpg(:) + frmean(:) - frame_accln(:))*dt !+ grav(:) 

		if(move_particles)then
			frame_vel(:) = frame_vel(:) + frame_accln(:)*dt
			frame_pos(:) = frame_pos(:) + frame_vel(:)*dt
		end if
    
		if (I_AM_NODE_ZERO) then
			WRITE(*,'(A25,3(2x,g17.8))')'FR MEAN    @ N  = ', (FRMEAN(j), j = 1, ndim)
			WRITE(*,'(A25,3(2x,g17.8))')'MPG        @ N  = ', (MPG(j), j = 1, ndim)
			if (move_particles) then
				WRITE(*,'(A25,3(2x,g17.8))')'FRAME ACCLN@ N  = ', (frame_accln(j), j = 1, ndim)
				WRITE(*,'(A25,3(2x,g17.8))')'FRAME VEL  @ N  = ', (frame_vel(j), j = 1, ndim)
				WRITE(*,'(A25,3(2x,g17.8))')'FRAME POS  @ N  = ', (frame_pos(j), j = 1, ndim)
			endif
			WRITE(*,'(A25,3(2x,g17.8))')'UMEAN      @ N+1= ', (UMEAN(j), j = 1, ndim)
		endif
	end subroutine velstep
end module steptheflow
