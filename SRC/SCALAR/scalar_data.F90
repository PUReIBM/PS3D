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


!GLOBAL DATA FOR SCALAR MODULES
!CURRENTLY SCALAR ROUTINES USE A LOT OF ARRAYS FROM THE HYDRODYNAMIC SIDE TO SAVE MEMORY.
!BE CAREFUL!

MODULE scalar_data
#include "../FLO/ibm.h"

  USE precision  
  USE constants 
  USE global_data,  ffphi=>ff
  IMPLICIT NONE 
  SAVE 
  
  INTEGER, PARAMETER :: nscalmax = 1 !maximum number of scalars possible
  REAL(PRCN), dimension(:), ALLOCATABLE :: sourcesink
  COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE ::  phif,nlphif,onlphif,sourcephif
  !COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE ::  ffphi
  COMPLEX(prcn), DIMENSION(:,:,:), ALLOCATABLE :: phiin,phiout
  !       Scalar real data 
  !REAL(prcn) :: phireal(mx, my, mz,nspmx)
  REAL(prcn), DIMENSION(:,:,:,:), ALLOCATABLE :: sourcephi
  !REAL(prcn), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gradphi
  !REAL(prcn) :: phi_an(mx,my,mz,nspmx)
  REAL(prcn), DIMENSION(:,:,:), ALLOCATABLE :: surf_scal_value, Nu3
  REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: phisurfall, phisurfall_nm1,phisurfall_nm2, &
          flux_body, flux_body2
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: phirmean,fphirmean,sum_flux_nm1,sum_flux_nm2,&
    sum_flux, GAMMA, phi_fluid_mean, phimodmean, phi_solid_mean
  
  REAL(PRCN), DIMENSION(:) , ALLOCATABLE :: nu_error_array, flux_global, flux_global2
  
  
  REAL(prcn) xis,xf,yi,yf,zi,zf,cons,uo,vo, nu_old, nu_error_hist, nu_error, phimean_des

  INTEGER ::  nxi,nxf,nyi,nyf,nzi,nzf,diff,nspmx, iter_scal
  REAL(prcn) :: phisurf,phistream,  Pr_or_Sc_nu!,gamma(nspmx)
  LOGICAL :: sourcepresent, LUMPED_CAP, zero_flow_in_solid, dorpr_scal, setphimean
  
  
  !Integer :: inn_itn_step

  NAMELIST/scal_propt/ nspmx,phistream, setphimean, phisurf,sourcepresent, LUMPED_CAP, &
    Pr_or_Sc_nu, zero_flow_in_solid, dorpr_scal

END MODULE scalar_data

