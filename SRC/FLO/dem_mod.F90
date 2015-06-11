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


MODULE dem_mod
  USE precision 
  USE constants 
  IMPLICIT NONE 
  SAVE
  PUBLIC
  DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

  INTEGER :: DIMN 
  !   Particle attributes
  !     Radius, density, mass, moment of inertia      
  INTEGER, PARAMETER:: DIM_M = 10
  INTEGER :: MMAX
  DOUBLE PRECISION S_TIME, DES_SPX_TIME, DES_RES_TIME, OVERLAP_MAX
  DOUBLE PRECISION DTSOLID, DTSOLID_FACTOR , lid_vel
  DOUBLE PRECISION P_TIME, PTC
  INTEGER NFACTOR, PARTICLES, NEIGH_MAX
  DOUBLE PRECISION AVG_RAD, RMS_RAD,  XLENGTH, YLENGTH, ZLENGTH, min_radius, max_radius, DES_F, DES_GAMMA, TSTOP
  LOGICAL, SAVE :: DES_ALLOC_CALLED=.FALSE., GENER_CONFIG_CASE = .FALSE.
  !
  !   Particle treatment at the walls  
  LOGICAL WALLFIXEDOVERLAP
  LOGICAL WALLDTSPLIT
  LOGICAL WALLREFLECT
  !
  !   Periodic Wall BC
  LOGICAL DES_PERIODIC_WALLS
  LOGICAL DES_PERIODIC_WALLS_X
  LOGICAL DES_PERIODIC_WALLS_Y
  LOGICAL DES_PERIODIC_WALLS_Z
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: IS_MOBILE, CAUSE_MOTION

  INTEGER :: NWALLS, TEST_PART

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RO_S, D_p0, test_force
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_RADIUS ! (PARTICLES)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RO_Sol ! (PARTICLES)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PVOL !(PARTICLES)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PMASS ! (PARTICLES)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMOI ! (PARTICLES)
  !
  !   Old and new particle positions, velocities (translational and
  !                                                             rotational)      
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_OLD ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_NEW ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_OLD ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_NEW ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_OLD ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_NEW ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPOS ! (PARTICLES,DIMN)
  !
  !   Total, normal and tangetial forces      
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FC ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FN ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FT ! (PARTICLES,DIMN)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FNS2 ! (DIMN)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FTS2 ! (DIMN)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FNS1 ! (DIMN)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FTS1 ! (DIMN)
  !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GRAV ! (DIMN)
  !
  !   Torque      
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TOW ! (PARTICLES,DIMN)
  !
  !   Accumulated spring forces      
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFN ! (PARTICLES,DIMN,MAXNEIGHBORS)
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFT ! (PARTICLES,DIMN,MAXNEIGHBORS)
  !
  !   Wall position, velocity and normal vector
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_POS ! (NWALLS,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_VEL ! (NWALLS,DIMN)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WALL_NORMAL ! (NWALLS,DIMN)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: PN ! (PARTICLES, MAXNEIGHBORS)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: PV ! (PARTICLES, MAXNEIGHBORS)
  !

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: NEIGHBOURS ! (PARTICLES, MAXNEIGHBORS)

  DOUBLE PRECISION RADIUS_EQ, NEIGHBOR_SEARCH_N
  DOUBLE PRECISION NEIGHBOR_SEARCH_RAD_RATIO, NEIGHBOR_SEARCH_DIST

  REAL(PRCN), DIMENSION(:,:,:,:), ALLOCATABLE  :: cgrid

  !   foctor for sum of radii in des_grid_based_neighbor_search
  DOUBLE PRECISION FACTOR_RLM
  !


  !   Particle-particle and Particle-wall contact parameters
  !             Spring contants      
  DOUBLE PRECISION KN, KN_W  ! Normal
  DOUBLE PRECISION KT, KT_W  ! Tangential
  !             Damping coeffients      
  DOUBLE PRECISION ETA_DES_N, ETA_N_W  ! Normal
  DOUBLE PRECISION ETA_DES_T, ETA_T_W  ! Tangential
  !             Damping coeffients in array form 
  DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE :: DES_ETAN, DES_ETAT !(MMAX, MMAX)

  DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE :: DES_ETAN_WALL, DES_ETAT_WALL !(MMAX)
  !             Friction coefiicients and coeff of restitution
  DOUBLE PRECISION MEW, MEW_W
  !coeff of restituion input in one D array, solid solid
  DOUBLE PRECISION DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2),DES_ET_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
  !coeff of restituion input in one D array, solid wall 
  DOUBLE PRECISION  DES_EN_WALL_INPUT(DIM_M),  DES_ET_WALL_INPUT(DIM_M)
  !actual coeff of rest.'s rearranged 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  REAL_EN, REAL_ET !(MMAX,MMAX)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  REAL_EN_WALL,  REAL_ET_WALL !(MMAX)
  !
  !   Wall treatment      
  INTEGER WALLCONTACT
  DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2
  !             Wall vibration parameters
  LOGICAL TEST_YMAXVAL, SHRINK
  DOUBLE PRECISION YMAXVAL


  Integer :: cx, cy,cz , MN
  REAL(prcn) :: dxc, dyc, dzc, particles_factor, cgrid_fac, dths, DTSOLID_ORIG
 
 
  INTEGER MAXNEIGHBORS


  !   Particles in a computational cell (for volume fraction)
  INTEGER, DIMENSION(:), ALLOCATABLE :: PINC 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIJK ! (PARTCILES,5)=>I,J,K,IJK,M 

  Integer, Dimension(:,:,:), Allocatable:: cnd
  REAL(prcn), Dimension(:), Allocatable:: XE, YN, ZT

  INTEGER :: IMIN1, IMAX1, JMIN1, JMAX1, KMIN1, KMAX1,  IMAX2, JMAX2, KMAX2, IMAX, JMAX, KMAX, IMAX3, JMAX3, KMAX3, IMIN2, JMIN2, KMIN2
  TYPE iap1
     INTEGER, DIMENSION(:), POINTER:: p
  END TYPE iap1
  !id's of particles in a cell 
  TYPE(iap1), DIMENSION(:,:,:), ALLOCATABLE:: pic

CONTAINS

  
  SUBROUTINE DES_CROSSPRDCT (AA, XX,YY) 
      IMPLICIT NONE
! 
      DOUBLE PRECISION AA(DIMN), XX(DIMN), YY(DIMN) 
! 
      IF(DIMN.EQ.3) THEN
         AA(1) = XX(2)*YY(3) - XX(3)*YY(2) 
         AA(2) = XX(3)*YY(1) - XX(1)*YY(3) 
         AA(3) = XX(1)*YY(2) - XX(2)*YY(1)
      ELSE
         AA(1) = - XX(1)*YY(2) 
         AA(2) = XX(1)*YY(1)  
      END IF

      RETURN  
      END SUBROUTINE DES_CROSSPRDCT


  

END MODULE dem_mod

  DOUBLE PRECISION FUNCTION DES_DOTPRDCT(XX,YY) 
    USE precision
    USE constants
    
    IMPLICIT NONE
    !
    
    INTEGER II
    INTEGER, PARAMETER :: DIMN=3
    DOUBLE PRECISION DOTP, XX(DIMN), YY(DIMN) 
    ! 
    DOTP = ZERO

    DO II = 1, DIMN
       DOTP = DOTP + XX(II)*YY(II)
    END DO
    DES_DOTPRDCT = DOTP
    
    RETURN  
  END FUNCTION DES_DOTPRDCT
