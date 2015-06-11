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


MODULE meanflo
#include "ibm.h"
  USE precision 
  USE constants 
  USE global_data
  USE tridiagonal
  USE gauss_seidel
  IMPLICIT NONE 
CONTAINS

  SUBROUTINE calc_meanflo(rks)
    
    !       this is a 3d code
    
    IMPLICIT NONE 
    INTEGER, INTENT(in) ::  rks
    
    !-----------------------------------------------------------------------
    !	local variables
    
    INTEGER :: j
    
    umean(1) = umean(1) + (-mpg(1)+frmean(1))*dt
    umean(2) = umean(2) + (-mpg(2)+frmean(2))*dt
    umean(3) = umean(3) + (-mpg(3)+frmean(3))*dt
    
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A25,3(2x,g12.5))')'FR MEAN   = ', (FRMEAN(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g12.5))')'MPG   = ', (MPG(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g12.5))')'UMEAN = ', (UMEAN(j), j = 1, ndim)
    end if
    RETURN
  end SUBROUTINE calc_meanflo
end MODULE meanflo


    
