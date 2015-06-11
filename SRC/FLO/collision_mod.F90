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


MODULE collision_mod
#include "ibm.h"
  USE precision 
  USE constants 
  USE global_data
  USE general_funcs
  USE soft_spring
  USE dem_mod

  IMPLICIT NONE 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MINCOLLS
  INTEGER, PRIVATE :: IFAC, FACTOR
Contains
  SUBROUTINE des_time_march(initialize)
    IMPLICIT NONE 
    LOGICAL, INTENT(in):: initialize
    INTEGER :: LL, NP, I, shrink_trial, max_shrink_trials, L
    logical :: test_ymaxval_tmp, shrink_tmp,free_evol_during_shrink
    REAL(prcn) :: ymin, ymax, Ay, By, tstop_tmp
    REAL(prcn):: ths, thcoll

    max_shrink_trials = 5
    shrink_trial = 0
    test_ymaxval_tmp = test_ymaxval
    tstop_tmp = tstop
    shrink_tmp = shrink
    free_evol_during_shrink = .true.
    if(I_AM_NODE_ZERO) CALL screen_separator(80,'C')
    if(I_AM_NODE_ZERO) WRITE(*,'(A20)')'IN COLLISION MODULES' 
    !    WRITE(*,*) 'IN DES_TIME MARCH, iniiti, generpartconfig =', initialize, gener_config_case

    IF(initialize)then

!!$  INITIALIZE COARSE GRID FOR DEM MODULES. The
!!$    particles will be moved based on this coarse grid.
       CALL INITIALIZE_COARSE_GRID

!!$Allocate memory for the arrays used in DEM module
       IF(.NOT.DES_ALLOC_CALLED) CALL DES_ALLOCATE_ARRAYS
!!$Initialize DEM arrays
       CALL DES_INIT_ARRAYS
       
!!$Assign Particle properties
       CALL CFASSIGN
       if(.not.GENER_CONFIG_CASE)then
          if(TRIM(collision_type).eq."softsphere")then
             DTSOLID_ORIG = DTSOLID
          end if
          if(I_AM_NODE_ZERO)then
             WRITE(*,'(A25,2(2x,g12.5))')'DES_EN_INPUT =', DES_EN_INPUT(1),DES_ET_INPUT(1)
             Write(*,'(2(A30,(2x,g17.8)))')'DT FLUID = ', DT, 'DT COLLISIONAL&
                  & = ', DTSOLID_ORIG
          end if
       end if

       CALL FIND_CELL_INDEX
       CALL PARTICLES_IN_CELL
       CALL GRID_BASED_NEIGHBOR_SEARCH
    END IF

    IF(.NOT.initialize)THEN
       IF(GENER_CONFIG_CASE) THEN
10000     continue          

          CALL  init_particles_jn
          IF(TRIM(collision_type).eq."eventdriven")THEN
             
             ALLOCATE( MINCOLLS(PARTICLES))
             MINCOLLS = 0
             do while (MINVAL(MINCOLLS(1:PARTICLES)).lt.min_colisions)
                CALL HARD_SPHERE_COLLISION
                CALL FIND_CELL_INDEX
                CALL PARTICLES_IN_CELL
                CALL GRID_BASED_NEIGHBOR_SEARCH
             end do
             DEALLOCATE(MINCOLLS)

          ELSE IF(TRIM(collision_type).eq."softsphere") THEN
             WRITE(*,*)'GENERATING INITIAL CONFIG BY SOFT SPHERE'
             WRITE(*,'(A25,2x,g12.5)')'DES_EN_INPUT =', DES_EN_INPUT(1)
             WRITE(*,*)'SHRINK, TESTYMAXVAL = ', SHRINK, TEST_YMAXVAL
             FACTOR = NINT(TSTOP/DTSOLID)
             WRITE(*,*)'FACTOR = ', FACTOR, TSTOP, DTSOLID
             S_TIME = ZERO 
             DO IFAC = 1, FACTOR
                !PRINT*,'SSCOLL, SHRINK, t =', SHRINK, S_TIME
                CALL CFUPDATEOLD

                CALL SOFT_SPHERE_COLLISION
                IF(TEST_YMAXVAL.AND.(.NOT.SHRINK)) EXIT
                CALL FIND_CELL_INDEX
                CALL PARTICLES_IN_CELL
                IF(MOD(IFAC,INT(NEIGHBOR_SEARCH_N)).EQ.0) CALL GRID_BASED_NEIGHBOR_SEARCH
             end DO

             IF(free_evol_during_shrink) then 
                test_ymaxval = test_ymaxval_tmp
                shrink = shrink_tmp
                tstop = tstop_tmp
             end IF


             IF(TEST_YMAXVAL.AND.SHRINK) THEN 


                IF(.not.free_evol_during_shrink)  then 

                   shrink_trial = shrink_trial + 1 

                   IF(shrink_trial.gt.max_shrink_trials) then 
                      WRITE(*,'(2(A20,2x,i4))') 'SHRINK TRIAL #', shrink_trial, 'GT THAN', MAX_SHRINK_TRIALS
                      WRITE(*,'(A)') 'FREELY EVOLVING THE SYSTEM FOR ONE STEP'

                      !Let the system evolve freely oncle
                      test_ymaxval = .false.
                      shrink = .false. 
                      shrink_trial = 0
                      free_evol_during_shrink = .true.


                      tstop = 3.d0*tstop 

                      ymin = MINVAL(DES_POS_NEW(1:PARTICLES,2)) - MAX_RADIUS
                      ymax = MAXVAL(DES_POS_NEW(1:PARTICLES,2)) + MAX_RADIUS

                      YLENGTH = ymax - ymin !+ 1.d0*MAX_RADIUS
                      Ay = YLENGTH/(ymax-ymin)
                      By = -Ay*ymin


                      DO L = 1, PARTICLES 
                         DES_POS_NEW(L,2) = Ay*DES_POS_NEW(L,2) + By
                      end DO

                      XC_GENER(1:PARTICLES,1:3) = DES_POS_NEW(1:PARTICLES, 1:3)
                      RAD_GENER(1:PARTICLES) = DES_RADIUS(1:PARTICLES)


                      DES_EN_INPUT(:) = 0.5
                      DES_EN_WALL_INPUT(:) = 1.0


!!$  INITIALIZE COARSE GRID FOR DEM MODULES. The
!!$    particles will be moved based on this coarse grid.
                      CALL INITIALIZE_COARSE_GRID

!!$Allocate memory for the arrays used in DEM module
                      IF(.NOT.DES_ALLOC_CALLED) CALL  DES_ALLOCATE_ARRAYS
!!$Initialize DEM arrays
                      CALL DES_INIT_ARRAYS

!!$Assign Particle properties
                      CALL CFASSIGN
                      CALL FIND_CELL_INDEX
                      CALL PARTICLES_IN_CELL
                      CALL GRID_BASED_NEIGHBOR_SEARCH


                      WRITE(*,'(A,2x,g17.8)') 'NEW YLENGTH = ', YLENGTH

                   ELSE


                      WRITE(*,'(A,i3,A)') 'SHRINKAGE NOT YET ACHEIVED AT',shrink_trial,' SO REDOING WITH NEW VELCOTIES'
                   end IF


                ELSE
                   WRITE(*,'(A,i3,A)') 'FINISHED THE FREELY EVOLVING STEP DURING SHRINK: NOW BACK TO SHRINK'
                   DES_EN_INPUT = 0.3
                   DES_EN_WALL_INPUT = 1.0


                   XC_GENER(1:PARTICLES,1:3) = DES_POS_NEW(1:PARTICLES, 1:3)
                   RAD_GENER(1:PARTICLES) = DES_RADIUS(1:PARTICLES)

!!$  INITIALIZE COARSE GRID FOR DEM MODULES. The
!!$    particles will be moved based on this coarse grid.
                   CALL INITIALIZE_COARSE_GRID

!!$Allocate memory for the arrays used in DEM module
                   IF(.NOT.DES_ALLOC_CALLED) CALL  DES_ALLOCATE_ARRAYS
!!$Initialize DEM arrays
                   CALL DES_INIT_ARRAYS

!!$Assign Particle properties
                   CALL CFASSIGN
                   CALL FIND_CELL_INDEX
                   CALL PARTICLES_IN_CELL
                   CALL GRID_BASED_NEIGHBOR_SEARCH


                   free_evol_during_shrink = .false.
                end IF


                goto 10000
             end IF

             IF(TEST_YMAXVAL) THEN 
                OPEN(1001, file=TRIM(RUN_NAME)//"_xc_post_ss_shrink.dat", form="formatted")

                write(1001,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" '
                do NP = 1, PARTICLES
                   WRITE(1001,'(10(2x,g15.8))')( DES_POS_NEW(NP, i), i = 1, dimn), DES_RADIUS(NP)
                ENDDO
                CLOSE(1001,status="keep")
             ELSE

                OPEN(1001, file=TRIM(RUN_NAME)//"_xc_post_ss_mfp.dat", form="formatted")
                write(1001,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" '
                do NP = 1, PARTICLES
                   WRITE(1001,'(10(2x,g15.8))')( DES_POS_NEW(NP, i), i = 1, dimn), DES_RADIUS(NP)
                ENDDO
                CLOSE(1001,status="keep")
             ENDIF
          end IF

       ELSE
          !DTSOLID = DT
          IF(collision_type.eq."eventdriven")THEN
             WRITE(*,'(A20)')'PERFORMING EVENTDRIVEN HARDSPHERE COLLISIONS'
             ths = dt
             thcoll = zero
             Write(*,'(A40,(2x,g17.6))')'CURRENT FLOW TIME STEP = ', dt
             Do While(thcoll.lt.ths)
                dths = ths - thcoll
                Write(*,'(A40,(2x,g17.6))')'REMAINING TIME IN HARD SPHERE = ', DTHS
                DTSOLID = DTHS
                CALL HARD_SPHERE_COLLISION
                CALL FIND_CELL_INDEX
                CALL PARTICLES_IN_CELL
                CALL GRID_BASED_NEIGHBOR_SEARCH
                CALL CFUPDATEOLD
                thcoll = thcoll + dtsolid
                WRITE(*,'(A40,2(2x,g17.6))')'CURRENT TIME IN HARD SPHERE = ', thcoll, ths-thcoll
             End Do

          ELSE IF(collision_type.eq."softsphere") THEN
             FACTOR = INT(DT/DTSOLID_ORIG)

				!JUST ADDED THIS TO LIMIT THE COLLISIONAL TIME STEP, ESPECIALLY FOR
				!ZERO SLIP CASES.
				if (FACTOR<50) then
					factor = 50
					dtsolid_ORIG = dt/factor
				elseif (FACTOR>200) then
					factor = 200
					dtsolid_ORIG = dt/factor
				endif
             
             DTSOLID = DTSOLID_ORIG
             if(FACTOR.LT.1)then 
                FACTOR = 1
                DTSOLID = DT
             end if

				if (I_AM_NODE_ZERO) then
					Write(*,'((A,3(2x,g12.5)))')'DT SOLID = ', DTSOLID
					WRITE(*,'((A,3(2x,I8)))')'NO. OF SS TIME STEPS = ', FACTOR
				endif             

             DO IFAC = 1, FACTOR
                test_force(:) = zero
                CALL SOFT_SPHERE_COLLISION
                CALL CFUPDATEOLD
                CALL FIND_CELL_INDEX
                CALL PARTICLES_IN_CELL
                !IF(MOD(IGLOBSTEP,INT(NEIGHBOR_SEARCH_N)).EQ.0) CALL GRID_BASED_NEIGHBOR_SEARCH
                CALL GRID_BASED_NEIGHBOR_SEARCH
                !WRITE(*,'(A,3(2x,g17.8)')'TOTAL FORCE = ', TEST_FORCE(:)
                !READ(*,*)
             END DO
          ELSE IF(collision_type.eq."none")THEN
             DTSOLID = dt
             DO LL = 1, PARTICLES 
                WRITE(*,'(2(A,3(2x,g12.5)))')'FORCE = ', RHOF*FORCE(LL,:)/PMASS(LL),' GRAV = ',(ONE-RHOF/RHOS)*GRAV(:)
                WRITE(*,'((A,3(2x,g12.5)))')'TORQUE = ', RHOF*OMOI(LL)*TORQ(LL,:)

! COMMENTING THIS, BECAUSE GRAVITY IS CHANGED TO MPG. IF BOTH GRAVITY AND MPG ARE AVAILABLE, THEN UNCOMMENT THIS
                DES_VEL_NEW(LL,:) = RHOF*FORCE(LL,:)/PMASS(LL) !+ (one-rhof/rhos)*GRAV(:) 

                DES_VEL_NEW(LL,:) = DES_VEL_OLD(LL,:) !+
                ! DES_VEL_NEW(LL,:)*DTSOLID
                OMEGA_NEW(LL,:) = OMEGA_OLD(LL,:) !+ RHOF*TORQ(LL,:)
                !*OMOI(LL)*DTSOLID
                
                Write(*,'(A30,3(2x,g17.8))')'DES_VEL_NEW =',&
                     & DES_VEL_NEW(LL,:)
                
                !                Write(*,'(A,3(2x,g17.8))')'DES_VEL_OLD =', DES_VEL_OLD(:)
                Write(*,'(A30,3(2x,g17.8))')'OMEGA_NEW =',&
                     & OMEGA_NEW(LL,:)
                CALL CFNEWVALUES(LL)
                CALL CFUPDATEOLD
                Write(*,'(A30,3(2x,g17.8))')'DES_POS_NEW =',&
                     & DES_POS_NEW(LL,:)/dx + one
                !READ(*,*)
             END DO
          END IF

!!$          CALL CFUPDATEOLD
!!$
!!$          CALL FIND_CELL_INDEX
!!$          CALL PARTICLES_IN_CELL
!!$          CALL GRID_BASED_NEIGHBOR_SEARCH
          
          CALL IBMUPDATE
       end IF !GENER_CONFIG_CASE

    END IF
    if(I_AM_NODE_ZERO)then
       Write(*,'(A20)')'LEAVING COLLISION MODULE'
       CALL screen_separator(80,'C')
    endif
!!$    DO I = 1, NBODY
!!$        PRINT*,'NEIGHS FOR I= ', I, ' ARE', NEIGHBOURS(I,1), NEIGHBOURS(I,2:NEIGHBOURS(I,1)+1)
!!$    END DO

  END SUBROUTINE des_time_march

  subroutine  initialize_coarse_grid
    implicit none 

    INTEGER ::  i,j,k, ii, ll

    REAL(prcn) :: dxeff

    MMAX = 1
    DIMN = NDIM
    NWALLS = 2*DIMN
    PARTICLES = NBODY
    NEIGHBOR_SEARCH_N = 10
    cgrid_fac = 1.5d0
    if(nphases.gt.1)then
       MN = 60
    else
       MN = 30
    end if
    FACTOR_RLM = 1.2d0
    PARTICLES_FACTOR = 1.0
    ! CX,CY, CZ are the number of nodes of the coarse grid excluding ghost nodes

    IF(GENER_CONFIG_CASE) THEN
       dxeff = two*(MAXVAL(RAD_GENER(1:NBODY)))*cgrid_fac
       !MAX_RADIUS = MAXVAL(RADBDY(1:NBODY))
    ELSE
       dxeff = two*(MAXVAL(RADBDY(1:NBODY)))*cgrid_fac*dx

       !MAX_RADIUS = MAXVAL(RADBDY(1:NBODY))*
    end IF

    cy = MAX(NINT(YLENGTH/dxeff),2)
    cx = MAX(NINT(XLENGTH/dxeff),2)
    cz = MAX(NINT(ZLENGTH/dxeff),2)

    !cz = cy
    !cx = MAX(NINT(doml(1)/doml(2))*cy,2)

    KN_W = 800000
    KN = 800000
    if(.not.GENER_CONFIG_CASE.and.(TRIM(collision_type).eq.'softsphere'))KN = 800000!20000

		!if (zero_slip) then
		!	if (ReT>small_number) then
		!		kn = rhos * pi/6 * (ReT*vis)**2 / dia_phys / (0.01)**2 ! the last term is the maximum overlap

		!		if (I_AM_NODE_ZERO) write (*,*) "KN CHANGED TO ", kn, " TO LIMIT MAXIMUM OVERLAP TO 0.01"
		!	endif
		!endif

    MEW = 0.0
    MEW_W = 0.0
    DTSOLID_FACTOR = 1.d0
    IF(allocated(cgrid))  then
       DEALLOCATE(cgrid)
    end IF

    ALLOCATE(cgrid(cx,cy,cz,3))

    dyc = YLENGTH/REAL(cy-1,prcn)
    dxc = XLENGTH/REAL(cx-1,prcn)
    dzc = ZLENGTH/REAL(cz-1,prcn)
    DO k = 1, cz
       do j = 1, cy
          do i = 1, cx 
             cgrid(i,j,k,1) = (i-1)*dxc
             cgrid(i,j,k,2) = (j-1)*dyc
             cgrid(i,j,k,3) = (k-1)*dyc
          end do
       end do
    end DO
    !PRINT*,'doml ', doml(1), doml(2)
    if (I_AM_NODE_ZERO) PRINT*,'cx, cy, cz = ', cx,cy,cz

    !PRINT*,'dxc, dyc, dzc = ', dxc,dyc,dzc
    ! IMAX,JMAX,KMAX are the number of physical cells, excluding the ghost cells
    IMAX = CX - 1
    JMAX = CY - 1
    KMAX = CZ - 1

    ! IMAX1,JMAX1,KMAX1 are the indices of the last physical CELL
    IMAX1 = IMAX + 1
    JMAX1 = JMAX + 1
    KMAX1 = KMAX + 1

    IMAX2 = IMAX1+1
    JMAX2 = JMAX1+1
    KMAX2 = KMAX1+1

    IMAX3 = IMAX2
    JMAX3 = JMAX2
    KMAX3 = KMAX2

    ! IMIN1,JMIN1,KMIN1 are the indices of the first physical CELL
    IMIN1 = 2
    JMIN1 = 2
    KMIN1 = 2

    IMIN2 = 1
    JMIN2 = 1
    KMIN2 = 1

!!$    PRINT*,'IMAX, JMAX, KMAX = ', IMAX, JMAX, KMAX
!!$    PRINT*,'IMIN1, JMIN1, KMIN1 = ', IMIN1, JMIN1, KMIN1
!!$    PRINT*,'IMIN2, JMIN2, KMIN2 = ', IMIN2, JMIN2, KMIN2
!!$
!!$    PRINT*,'IMAX1, JMAX1, KMAX1 = ', IMAX1, JMAX1, KMAX1
!!$    PRINT*,'IMAX2, JMAX2, KMAX2 = ', IMAX2, JMAX2, KMAX2
!!$    PRINT*,'IMAX3, JMAX3, KMAX3 = ', IMAX3, JMAX3, KMAX3

    IF(ALLOCATED(XE)) THEN

       DEALLOCATE(XE,YN, ZT)
       DEALLOCATE(pic, cnd)
    end IF

    ALLOCATE(XE(IMAX2), YN(JMAX2), ZT(KMAX2))

    ALLOCATE(pic(IMAX2,JMAX2,KMAX2))
    ALLOCATE(cnd(IMAX2,JMAX2,KMAX2))

    DO  k  = 1,KMAX2!MAX(KMAX1-1,1)
       DO j  = 1,JMAX2
          DO  i  = 1,IMAX2
             NULLIFY(pic(i,j,k)%p)
             cnd(i,j,k) = 0
          end DO
       end DO
    end DO


    XE = ZERO
    YN = ZERO
    ZT = ZERO
    XE(1) = ZERO
    YN(1) = ZERO
    ZT(1) = ZERO
    DO I = IMIN1, IMAX2
       XE(I) = XE(I-1) + DXC
       !PRINT*,'XE = ', I, XE(I)
    END DO

    DO J  = JMIN1, JMAX2
       YN(J) = YN(J-1) + DYC
    END DO

    DO K = KMIN1, KMAX2
       ZT(K) = ZT(K-1) + DZC
    END DO
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A40,2(2x,g17.8))') 'XLENGTH IN AND CALC =', XLENGTH, XE(IMAX1) - XE(IMIN2)
       WRITE(*,'(A40,2(2x,g17.8))') 'YLENGTH IN AND CALC =', YLENGTH, YN(JMAX1) - YN(JMIN2)
       WRITE(*,'(A40,2(2x,g17.8))') 'ZLENGTH IN AND CALC =',  ZLENGTH, ZT(KMAX1) - ZT(KMIN2)
    end if
  end subroutine initialize_coarse_grid

  SUBROUTINE DES_ALLOCATE_ARRAYS 

    IMPLICIT NONE

    INTEGER NPARTICLES, I,J,K
    INTEGER :: DIMENSION_I, DIMENSION_J, DIMENSION_K 

    DES_ALLOC_CALLED = .TRUE.
    DIMENSION_I   = IMAX3
    DIMENSION_J   = JMAX3
    DIMENSION_K   = KMAX3
    
    !particles = npc*(imax)*(jmax)*kmax
    NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS
    MAXNEIGHBORS = MN + 1 + NWALLS
    ALLOCATE(RO_S(MMAX), D_p0(MMAX))

    Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )

    ALLOCATE(PIJK(PARTICLES,3))
    
    ALLOCATE(IS_MOBILE(PARTICLES), CAUSE_MOTION(PARTICLES))
    
    ALLOCATE(REAL_EN(MMAX,MMAX),REAL_ET(MMAX,MMAX))
    ALLOCATE(REAL_EN_WALL(MMAX),REAL_ET_WALL(MMAX))
    !
    ALLOCATE(DES_ETAN(MMAX,MMAX))
    ALLOCATE(DES_ETAT(MMAX,MMAX))
    ALLOCATE(DES_ETAN_WALL(MMAX), DES_ETAT_WALL(MMAX))
    Allocate(  DES_RADIUS (NPARTICLES) )
    Allocate(  RO_Sol (NPARTICLES) )
    Allocate(  PVOL (NPARTICLES) )
    Allocate(  PMASS (NPARTICLES) )
    Allocate(  OMOI (NPARTICLES) )
    !
    !   Old and new particle positions, velocities (translational and
    !                                                             rotational) )      
    Allocate(  DES_POS_OLD (NPARTICLES,DIMN) )
    Allocate(  DES_POS_NEW (NPARTICLES,DIMN) )
    Allocate(  DES_VEL_OLD (NPARTICLES,DIMN) )
    Allocate(  DES_VEL_NEW (NPARTICLES,DIMN) )
    IF(DIMN.GT.2) THEN
       Allocate(  OMEGA_OLD (NPARTICLES,DIMN) )
       Allocate(  OMEGA_NEW (NPARTICLES,DIMN) )
    ELSE
       Allocate(  OMEGA_OLD (NPARTICLES,1) )
       Allocate(  OMEGA_NEW (NPARTICLES,1) )
    END IF
    Allocate(  PPOS (NPARTICLES,DIMN) )
    !
    !   Total, normal and tangetial forces      
    Allocate(  FC (NPARTICLES,DIMN) )
    Allocate(  FN (NPARTICLES,DIMN) )
    Allocate(  FT (NPARTICLES,DIMN) )
    Allocate(  FNS2 (DIMN) )
    Allocate(  FTS2 (DIMN) )
    Allocate(  FNS1 (DIMN) )
    Allocate(  FTS1 (DIMN) )

    !
    !   Torque     
    IF(DIMN.EQ.3) THEN 
       Allocate(  TOW (NPARTICLES,DIMN) )
    ELSE
       Allocate(  TOW (NPARTICLES,1) )
    END IF
    !
    !   Accumulated spring forces      
    Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )
    Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )
    !
    !   Wall position, velocity and normal vector
    Allocate(  DES_WALL_POS (NWALLS,DIMN) )
    Allocate(  DES_WALL_VEL (NWALLS,DIMN) )
    Allocate(  WALL_NORMAL (NWALLS,DIMN) )
    Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
    Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )

    ALLOCATE(TEST_FORCE(DIMN))
    !   Particles in a computational cell (for volume fraction) )
    !Allocate(  PIJK (PARTICLES,5) )


  end SUBROUTINE DES_ALLOCATE_ARRAYS




  SUBROUTINE CFASSIGN

    IMPLICIT NONE
    LOGICAL:: filexist, isopen

    INTEGER L, IJK, M, I,J, K, COUNT_E
    DOUBLE PRECISION FOUR_BY_THREE, RAD2, MINMASS, MASS_I, MASS_J, MASS_EFF
    DOUBLE PRECISION :: TCOLL, TCOLL_TMP, AVG_MASS, MAXMASS

    !     
    !---------------------------------------------------------------------
    !     Assignments
    !---------------------------------------------------------------------
    !     
    WALLDTSPLIT = .FALSE.
    FOUR_BY_THREE = 4.0d0/3.0d0
    MINMASS = LARGE_NUMBER
    MAXMASS = SMALL_NUMBER
    MAX_RADIUS = ZERO
    MIN_RADIUS = LARGE_NUMBER
    TCOLL = LARGE_NUMBER
    RMS_RAD = ZERO
    !if(.not.gener_config_case)KN = 2.0E+04
    DO L = 1, PARTICLES

       if(GENER_CONFIG_CASE)then
          Ro_Sol(L) = three/(four*pi*DES_RADIUS(L)**3.d0)
       else
          Ro_Sol(L) = RHOS!three/(four*pi*DES_RADIUS(L)**3.d0)
       end if

       RAD2 = DES_RADIUS(L)**2
       PVOL(L) = FOUR_BY_THREE*Pi*RAD2*DES_RADIUS(L)
       PMASS(L) = PVOL(L)*RO_Sol(L) 
       OMOI(L) = 2.5d0/(PMASS(L)*RAD2) !one over MOI
       MAX_RADIUS = MAX(MAX_RADIUS, DES_RADIUS(L))
       MIN_RADIUS = MIN(MIN_RADIUS, DES_RADIUS(L))
       RMS_RAD = RMS_RAD + DES_RADIUS(L)**2.d0
       IF(PMASS(L).LT.MINMASS) MINMASS = PMASS(L) 
       MAXMASS = MAX(PMASS(L), MAXMASS)
       
    END DO
    RMS_RAD = SQRT(RMS_RAD/PARTICLES)
    AVG_MASS = SUM(PMASS(1:PARTICLES))/PARTICLES
    AVG_RAD = SUM(DES_RADIUS(1:PARTICLES))/PARTICLES
    KT = (2.d0/7.d0)*KN
    KT_W = (2.d0/7.d0)*KN_W
    
    IF(.NOT.XPERIODIC(1)) THEN
       DES_PERIODIC_WALLS = .FALSE.
       DES_PERIODIC_WALLS_X = .FALSE.
       DES_PERIODIC_WALLS_Y = .FALSE.
       DES_PERIODIC_WALLS_Z = .FALSE.
    ELSE
       DES_PERIODIC_WALLS = .TRUE.
       DES_PERIODIC_WALLS_X = .TRUE.
       DES_PERIODIC_WALLS_Y = .TRUE.
       DES_PERIODIC_WALLS_Z = .TRUE.
    ENDIF

    RADIUS_EQ = DES_RADIUS(1)*1.05D0
    NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO * RADIUS_EQ
    !DTSOLID = pi*SQRT(one/(KN/PMASS(1) - (ETA_DES_N**2)/4.d0))
    !DTSOLID = DTSOLID/50

    DTSOLID = DTSOLID_FACTOR*2.0D0*PI*SQRT((MINMASS)/(15*KN)) ! DTs - Rotational Constraint
#if 0
    if(.not.GENER_CONFIG_CASE)then
       PRINT*,' In CFASSIGN'
       PRINT*,'MINMASS = ', MINMASS
       PRINT*,'Kn = ', Kn
       PRINT*,'DTSOLID FACTOR = ', DTSOLID_FACTOR
       READ(*,*)
    end if
    !     DTSOLID = DTSOLID_FACTOR*2D0*PI*SQRT(MINMASS/(6*KN)) ! DTs - Translational Constraint

    !Print*,'DTSOLID = ', dtsolid
    !Print*,'MAX_RADIUS = ', MAX_RADIUS

    !read(*,*)
#endif
    WX1 = ZERO 
    EX2 = XLENGTH 
    BY1 = ZERO
    TY2 = YLENGTH 

    SZ1 = ZERO 
    NZ2 = ZLENGTH

    
    !ARRANGE THE COEFF OF RESTITUTION MATRIX FROM INPUT EN VALUES
    count_e = 0
    DO I = 1, MMAX
       DO J = I, MMAX
          COUNT_E = COUNT_E + 1
          
          
          REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
          REAL_ET(I,J) = DES_ET_INPUT(COUNT_E)
          
          MASS_I = (PI*(D_P0(I)**3.d0)*RO_S(I))/6.d0
          MASS_J = (PI*(D_P0(J)**3.d0)*RO_S(J))/6.d0
          MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)
          DES_ETAN(I,J) = 2.D0*SQRT(KN*MASS_EFF)*ABS(LOG(REAL_EN(I,J)))
          !PRINT*,'MASSI, MASSJ = ', MASS_I, MASS_J, MASS_EFF, KN

          DES_ETAN(I,J) = DES_ETAN(I,J)/SQRT(PI*PI + (LOG(REAL_EN(I,J)))**2.0)
          DES_ETAT(I,J) = HALF*DES_ETAN(I,J)

          !new TCOLL_TMP = PI/SQRT(KN/MASS_EFF - ((DES_ETAN(I,J)/MASS_EFF)**2.d0)/4.d0)
	  TCOLL_TMP = 1./SQRT(KN/MASS_EFF)
          !WRITE(*,*) 'KN, MASS EFF = ', KN, MASS_EFF, DES_ETAN(I,J)
          TCOLL = MIN(TCOLL_TMP, TCOLL)
       ENDDO
    ENDDO

    COUNT_E = 0 
    DO I = 1, MMAX
       COUNT_E = COUNT_E + 1  
       REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
       REAL_ET_WALL(I) = DES_ET_WALL_INPUT(COUNT_E)
       MASS_I = (PI*(D_P0(I)**3.d0)*RO_S(I))/6.d0
       MASS_J = MASS_I
       MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)

       DES_ETAN_WALL(I) = 2.d0*SQRT(KN_W*MASS_EFF)*ABS(LOG(REAL_EN_WALL(I)))
       DES_ETAN_WALL(I) = DES_ETAN_WALL(I)/SQRT(PI*PI + (LOG(REAL_EN_WALL(I)))**2.0)
       DES_ETAT_WALL(I) = HALF*DES_ETAN_WALL(I)

    ENDDO

    DO I = 1, MMAX
       DO J = I, MMAX
          REAL_EN(J, I) = REAL_EN(I,J)
          REAL_ET(J, I) = REAL_ET(I,J)
          DES_ETAN(J,I) = DES_ETAN(I,J)
          DES_ETAT(J,I) = DES_ETAT(I,J)
       ENDDO
    ENDDO

!!$    DO I = 1, MMAX
!!$       DO J = 1, MMAX
!!$          WRITE(*,*) 'I AND J = ', I, J
!!$          WRITE(*,*) 'REAL_EN AND ET  = ', REAL_EN(I,J), REAL_ET(I,J)
!!$       ENDDO
!!$    ENDDO
!!$
!!$
!!$    DO I = 1, MMAX
!!$       DO J = 1, MMAX
!!$          WRITE(*,*) 'I AND J = ', I, J
!!$          WRITE(*,*) 'ETA_N AND ETA_T  = ', DES_ETAN(I,J), DES_ETAT(I,J)
!!$       ENDDO
!!$    ENDDO

    DTSOLID = TCOLL/50.d0
    
    if(I_AM_NODE_ZERO)WRITE(*,*) 'MIN TCOLL AND DTSOLID = ', TCOLL, DTSOLID
    !READ(*,*)
    RETURN
  END SUBROUTINE CFASSIGN

  subroutine find_cell_index
    implicit none 
    integer :: ip
    real(prcn) :: tempx, tempy, tempz

    CND(:,:,:) = 0
    DO IP = 1, PARTICLES
       
       tempx = des_pos_new(IP,1)
       tempy = des_pos_new(IP,2)
       tempz = des_pos_new(IP,3)
       !PRINT*,'tempx = ', tempx, dxc,tempx/dxc, int(tempx/dxc)

       PIJK(IP,1) = MIN(INT(tempx/dxc)+2, IMAX1)
       PIJK(IP,2) = MIN(INT(tempy/dyc)+2, JMAX1)
       PIJK(IP,3) = MIN(INT(tempz/dzc)+2, KMAX1)
       !PRINT*,'PC: ', PIJK(IP,1), PIJK(IP,2), PIJK(IP,3)
       
       cnd(PIJK(IP,1),PIJK(IP,2),PIJK(IP,3)) = cnd(PIJK(IP,1),PIJK(IP,2),PIJK(IP,3)) + 1
       
    end DO
  end subroutine find_cell_index
  
  subroutine particles_in_cell
    implicit none 
    INTEGER :: PC(3), npic , i,j,k, ip, pos
    
    INTEGER, DIMENSION(1:cx,1:cy,1:cz):: icount
    
    DO  k = 1,KMAX2        !MAX(KMAX1-1,1)
       DO  j = 1,JMAX2
          DO  i = 1,IMAX2
             NPIC = CND(i,j,k)
             !PRINT*,'NPIC = ', NPIC
             IF (ASSOCIATED(pic(i,j,k)%p)) THEN
                IF (npic.NE.SIZE(pic(i,j,k)%p)) THEN
                   DEALLOCATE(pic(i,j,k)%p)
                   IF (npic.GT.0) then 
                      !PRINT*,'NPIC = ', NPIC, i,j,k
                      ALLOCATE(pic(i,j,k)%p(npic))
                   end IF
                   
                ENDIF
             ELSE
                IF(npic.GT.0) ALLOCATE(pic(i,j,k)%p(npic))
                             
             ENDIF
             
          end DO
       end DO
    end DO
    !    PRINT*,'CND = ', CND
    icount(:,:,:) = 1
    
    DO ip = 1, NBODY
       PC(:) =  PIJK(IP,1:3)
       !PRINT*,'PC = ', PC(:), SIZE(pic(PC(1), PC(2), PC(3))%p), CND(PC(1), PC(2), PC(3))
       pos = icount(pc(1),pc(2),pc(3))
       pic(pc(1),pc(2),pc(3))%p(pos) = ip
       icount(pc(1),pc(2),pc(3)) = &
            & icount(pc(1),pc(2),pc(3)) + 1
    ENDDO
  end subroutine particles_in_cell

  SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
    USE general_funcs
    USE global_data
    IMPLICIT NONE
    LOGICAL PER_COND, ALREADY_NEIGHBOURS
    INTEGER I, II, LL, CO, NI, TEMP, JJ, KK , J, K, NEIGH_L, L_MAX, PNO_MAX
    INTEGER KM1, KP1, IM1, IP1, JM1, JP1, PNO, NPG, PC(3), IP2, NLIM
    DOUBLE PRECISION  DIST(NDIM), DISTMAG, R_LM, LX, LY, LZ, XPER_FAC, YPER_FAC, ZPER_FAC,  CORD_PNO(NDIM), TMP_OVERLAP

    !DOUBLE PRECISION :: DES_DOTPRDCT 
    
    !CALL screen_separator(80,'-')
    !PRINT*, 'IN CELL LINKED LIST SEARCH'
    DO I = 1, NBODY
       DO II = 1, MAXNEIGHBORS
          NEIGHBOURS(I,II) = -1
       END DO
       NEIGHBOURS(I,1) = 0
    END DO


    LX = XE(IMAX1) - XE(1)
    LY = YN(JMAX1) - YN(1)
    LZ = ZT(KMAX1) - ZT(1)
    OVERLAP_MAX = SMALL_NUMBER
    PNO_MAX = 0
    L_MAX= 0
!!$    DO LL = 1, PARTICLES 
!!$       DES_POS_NEW(LL,1) = (XC(LL,1)+foffset-one)*dx
!!$       DES_POS_NEW(LL,2) = (XC(LL,2)-one)*dy
!!$       DES_POS_NEW(LL,3) = (XC(LL,3)-one)*dz
!!$       DES_RADIUS(LL) = RADBDY(LL)*dx
!!$    end DO

    !WRITE(*,'(A,/,3(2x,g12.5,/))') ' DOMAIN LENGTH = ', LX, LY, LZ
    DO LL = 1, PARTICLES
       NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
       NLIM  = NEIGHBOURS(LL,1) + 1
       IF(NLIM.GT.MAXNEIGHBORS) THEN 
          WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE ', LL 
          STOP 
       end IF
       
       NEIGHBOURS(LL,NLIM) = LL

       PC(:) =  PIJK(LL,1:3)!+1
       
       II = PC(1)
       JJ = PC(2)
       KK = PC(3)
       IP1 = II+1
       IM1 = II-1
       JP1 = JJ+1
       JM1 = JJ-1
       
       IF(DIMN.EQ.3) THEN 
          KP1 = KK+1
          KM1 = KK-1
       end IF
       
       
       DO KK = KM1, KP1
          DO JJ = JM1, JP1
             DO II = IM1, IP1
                
                I = II
                J = JJ
                K = KK
                
                XPER_FAC = 0
                YPER_FAC = 0
                ZPER_FAC = 0
                PER_COND = .FALSE.
                IF(II.GT.IMAX1) THEN 
                   IF(INTX_PER) THEN 
                      I = IMIN1
                      XPER_FAC = one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true EAST',I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      I = IMAX1
                   ENDIF
                ENDIF
                
                IF(II.LT.IMIN1) THEN 
                   IF(INTX_PER) THEN 
                      I = IMAX1
                      XPER_FAC = -one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true WEST', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE 
                      I = IMIN1
                   ENDIF
                ENDIF

                IF(JJ.GT.JMAX1) THEN 
                   IF(INTY_PER) THEN 
                      J = JMIN1
                      YPER_FAC = one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true NORTH', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      J = JMAX1
                   ENDIF
                ENDIF

                IF(JJ.LT.JMIN1) THEN 
                   IF(INTY_PER) THEN 
                      J = JMAX1
                      YPER_FAC = -one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true SOUTH', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      J = JMIN1
                   ENDIF
                ENDIF

                IF(DIMN.EQ.3) THEN 
                   IF(KK.GT.KMAX1) THEN 
                      IF(INTZ_PER) THEN
                         K = KMIN1
                         ZPER_FAC = one
                      ELSE
                         K = KMAX1
                      ENDIF
                   ENDIF

                   IF(KK.LT.KMIN1) THEN 
                      IF(INTZ_PER) THEN
                         K = KMAX1
                         ZPER_FAC = -one
                      ELSE
                         K = KMIN1
                      ENDIF
                   ENDIF
                ENDIF

                If (ASSOCIATED(PIC(I,J,K)%p)) then
                   NPG = SIZE(PIC(I,J,K)%p)
                Else
                   NPG = 0
                Endif


                Do IP2 = 1,NPG
                   PNO = PIC(I,J,K)%p(ip2)


                   if(PNO.GT.LL) then 

                      R_LM = DES_RADIUS(LL) + DES_RADIUS(PNO)!+1.5*dxc
                      R_LM = FACTOR_RLM*R_LM
                      CORD_PNO(1) = DES_POS_NEW(PNO,1) + XPER_FAC*(LX)
                      CORD_PNO(2) = DES_POS_NEW(PNO,2) + YPER_FAC*(LY)
                      IF(DIMN.EQ.3) THEN 
                         CORD_PNO(3) = DES_POS_NEW(PNO,3) + ZPER_FAC*(LZ)
                      ENDIF


                      DIST(:) = CORD_PNO(:) - DES_POS_NEW(LL,:)
                      DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))

                      ALREADY_NEIGHBOURS = .FALSE.

!!$                      IF(LL.EQ.1.AND.PNO.EQ.2) THEN 
!!$                         WRITE(*,*)'CORD=', CORD_PNO(1)
!!$                         Write(*,*)'POS1 : ', DES_POS_NEW(LL,2),DES_RADIUS(LL)
!!$                         Write(*,*)'POS2 : ', DES_POS_NEW(PNO,2), DES_RADIUS(PNO)
!!$                         WRITE(*,*)'DISTMAG', DISTMAG, R_LM!-DISTMAG
!!$                      ENDIF

                      DO NEIGH_L = 2, NEIGHBOURS(LL,1)+1
                         IF(PNO.EQ. NEIGHBOURS(LL,NEIGH_L)) ALREADY_NEIGHBOURS=.true.
                      ENDDO

                      IF(R_LM - DISTMAG.gt.SMALL_NUMBER.AND.(.NOT.ALREADY_NEIGHBOURS)) THEN 
                         TMP_OVERLAP = ((DES_RADIUS(LL) + DES_RADIUS(PNO))-DISTMAG)/(DES_RADIUS(LL) + DES_RADIUS(PNO))
                         TMP_OVERLAP = TMP_OVERLAP*100
                         IF(TMP_OVERLAP.GT.OVERLAP_MAX) THEN 

                            OVERLAP_MAX = MAX(OVERLAP_MAX, TMP_OVERLAP)
                            L_MAX = LL
                            PNO_MAX = PNO

                         end IF



!!$                           IF(PER_COND) THEN 
!!$                              WRITE(*,*) 'pC = ', pc
!!$                              WRITE(*,*) 'II, JJ = ', II, JJ
!!$                              WRITE(*,*) 'I, J = ', I, J
!!$                              WRITE(*,*) 'XYPER_FAC ', XPER_FAC, YPER_FAC
!!$                              WRITE(*,*) 'DES_VEL_NEW = ', DES_POS_NEW(PNO,:)
!!$                              WRITE(*,*) 'MODIFIED POSITION = ', CORD_PNO(:)
!!$                           ENDIF


                         NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
                         NLIM  = NEIGHBOURS(LL,1) + 1
                         IF(NLIM.GT.MAXNEIGHBORS) THEN 
                            if (I_AM_NODE_ZERO) WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE LL', LL 
                            if (I_AM_NODE_ZERO) WRITE(*,*) 'EITHER REDUCE THE R_LM FACTOR OR INCREASE MN IN MFIX.DAT'
                            if (I_AM_NODE_ZERO) PRINT*,'POSL = ',DES_POS_NEW(LL,:)
                            DO NEIGH_L = 2, NEIGHBOURS(LL,1)+1

                               DIST(:) = DES_POS_NEW(NEIGHBOURS(LL,NEIGH_L),:) - DES_POS_NEW(LL,:)
                               DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
                               if (I_AM_NODE_ZERO) PRINT*,'LL =',NEIGHBOURS(LL,NEIGH_L), DES_POS_NEW(NEIGHBOURS(LL,NEIGH_L),:)
                               if (I_AM_NODE_ZERO) PRINT*,DISTMAG, FACTOR_RLM*(DES_RADIUS(LL) + DES_RADIUS(NEIGHBOURS(LL,NEIGH_L))),DES_RADIUS(LL),  DES_RADIUS(NEIGHBOURS(LL,NEIGH_L)), FACTOR_RLM
                            ENDDO
                            STOP 
                         end IF
                         NEIGHBOURS(LL,NLIM) = PNO

                         NEIGHBOURS(PNO,1) = NEIGHBOURS(PNO,1) + 1
                         NLIM  = NEIGHBOURS(PNO,1) + 1
                         IF(NLIM.GT.MAXNEIGHBORS) THEN 
                            WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE PNO', PNO
                            WRITE(*,*) 'EITHER REDUCE THE R_LM FACTOR OR INCREASE MN IN MFIX.DAT'
                            STOP 
                         end IF
                         NEIGHBOURS(PNO,NLIM) = LL

                      end IF !contact condition
                   end if !PNO.GT.LL
                end Do !IP2
             end DO
          end DO
       end DO
    end DO

!!$    do LL = 1, PARTICLES 
!!$       PRINT*,'NEIGHBORS = ', NEIGHBOURS(LL,:)
!!$    end do


    IF(L_MAX.NE.0.AND.(.not.gener_config_case))THEN 
       if (I_AM_NODE_ZERO) then
			!CALL screen_separator(80,'@')
			PRINT*,'MAXIMUM OVERLAP = @', OVERLAP_MAX, L_MAX, PNO_MAX, '@'
			!CALL screen_separator(80,'@')       
			!Write(unit_overlap,'(I8,2x,g17.6)') totcolls, overlap_max
		endif


       DIST(:) = DES_POS_NEW(L_MAX, :) - DES_POS_NEW(PNO_MAX,:)
       DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
       
       !PRINT*,  L_MAX, PNO_MAX, DISTMAG, DES_RADIUS(L_MAX)+DES_RADIUS(PNO_MAX)

       !PRINT*,'MAXIMUM OVERLAP, PART POSs L', DES_POS_NEW(L_MAX,:)
       !PRINT*,'MAXIMUM OVERLAP, PART POSs J', DES_POS_NEW(PNO_MAX,:)
       !PRINT*,'MAXIMUM OVERLAP, PART CELLS L', PIJK(L_MAX,:)
       !PRINT*,'MAXIMUM OVERLAP, PART CELLS J', PIJK(PNO_MAX,:)
       !READ(*,*)

    end IF

  END SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
  

  SUBROUTINE compute_vol_frac_cgrid
    USE interpolation
    IMPLICIT NONE
    INTEGER :: i,j,k, n, m,l, pc(3), onew, ii,jj,kk &
         & , ib,ie,jb,je,kb,ke, im, jm, km

    REAL(prcn) ::  pl,nll(ndim),onll(ndim),ppll(ndim),dfll(ndim), pos(3), ul(3), xl(3)
    maxvolfrac = SMALL_NUMBER
    !weightbar = 0.d0

    im = cx
    jm = cy
    km = cz
    if (I_AM_NODE_ZERO) PRINT*,'In compue vol_fracr, interp_scheme = ', interp_scheme

    DO m=1,nbody          
       xl(1)=(xc(m,1))*dx
       xl(2)=(xc(m,2))*dy
       xl(3)=(xc(m,3))*dz
       pc(:) = PIJK(m,:)
       !WRITE(*,*)'PC = ',xl
       SELECT CASE(interp_scheme)

       CASE('lpi')
          !print*, 'order in set stencil = ', order
          ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
          ie = MIN(im,pc(1) + ob2r)
          if(.not.intx_per) then 
             IF (ib.EQ.1 ) ie = ib + order - 1
             IF (ie.EQ.im) ib = ie - order + 1
          else 
             IF (ib.EQ.1 ) ib = ie - order + 1
             IF (ie.EQ.im) ie = ib + order - 1
          end if

          jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
          je = MIN(jm,pc(2) + ob2r)
          if(.not.inty_per) then
             IF (jb.EQ.1 ) je = jb + order - 1
             IF (je.EQ.jm) jb = je - order + 1
          else
             IF (jb.EQ.1 ) jb = je - order + 1
             IF (je.EQ.jm) je = jb + order - 1
          end if

          kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
          ke = MIN(km,pc(3) + ob2r)
          If(.not.intz_per) then 
             IF (kb.EQ.1 ) ke = kb + order - 1
             IF (ke.EQ.km) kb = ke - order + 1
          else
             IF (kb.EQ.1 ) kb = ke - order + 1
             IF (ke.EQ.km) ke = kb + order - 1
          end If
       end SELECT

       onew = order

       do k = 1, onew
          do j = 1, onew
             do i = 1, onew
                ii = ib+i-1
                jj = jb+j-1
                kk = kb+k-1
                gstencil(i,j,k,1) = (ii-1)*dxc
                gstencil(i,j,k,2) = (jj-1)*dyc   
                gstencil(i,j,k,3) = (kk-1)*dzc
                if(ii.lt.1.and.intx_per) ii = cx+ii !new-1
                if(ii.gt.cx.and.intx_per) ii = ii-cx !new+1
                if(jj.lt.1.and.inty_per) jj = cy+jj
                if(jj.gt.cy.and.inty_per) jj = jj-cy
                if(kk.lt.1.and.intz_per) kk = cz+kk
                if(kk.gt.cz.and.intz_per) kk = kk-cz
!!$                  
                vsten(i,j,k,1:ndim) = zero

             end do
          end do
       end do
       CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:3),&
            & vsten(1:onew,1:onew,1:onew,1:ndim),xl(1:ndim),ul(1:ndim),onew,&
            & interp_scheme,weightp) 


       do k = 1, onew 
          do j = 1, onew
             do i = 1, onew
                DO n=1,ndim
                   ii = ib+i-1
                   jj = jb+j-1
                   kk = kb+k-1

                   if(ii.lt.1.and.intx_per) ii = cx+ii !new -1
                   if(ii.gt.cx.and.intx_per) ii = ii-cx !new +1
                   if(jj.lt.1.and.inty_per) jj = cy+jj
                   if(jj.gt.cy.and.inty_per) jj = jj-cy
                   if(kk.lt.1.and.intz_per) kk = cz+kk
                   if(kk.gt.cz.and.intz_per) kk = kk-cz
                   !weightbar(ii,jj,kk) = weightbar(ii,jj,kk) + (4.d0*pi*((radbdy(m)*dx)**3.d0)/(3.d0))*weightp(i,j,k)

                ENDDO
             ENDDO
          ENDDO
       ENDDO


    end DO!CLOSE LOOP OVER ALL BODIES

!!$
    do k = 1, cz
       do j = 1, cy
          do i = 1, cx
             !weightbar(i,j,k) = weightbar(i,j,k)/(dxc*dyc*dzc)
             !maxvolfrac = MAX(maxvolfrac, weightbar(i,j,k))
          end do
       end do
    end do
    if (I_AM_NODE_ZERO) PRINT*,'MAXVOLFRAC = ', maxvolfrac
    open(1000, file="volfrac.dat", form="formatted")

    write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "volfrac" '

    write(1000,*)'ZONE F=POINT, I=', cx,  ', J=', cy, ', K=', cz


!!$    do k = 1, cz
!!$       do j = 1, cy
!!$          do i = 1, cx
!!$                write(1000,'(3(2x,i4),2x,g12.5)')(i),(j),(k),weightbar(i,j,k)
!!$             enddo
!!$          enddo
!!$       enddo

    !       close(1000, status="keep")

  end SUBROUTINE compute_vol_frac_cgrid

  SUBROUTINE DES_INIT_ARRAYS

    IMPLICIT NONE
    !-----------------------------------------------
    !     G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !     L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !     L o c a l   V a r i a b l e s
    !-----------------------------------------------
    !     loop counters
    INTEGER :: M, N, K, LL
    !     
    !     Coefficient of restitution (old symbol)
    DOUBLE PRECISION :: E, TEST_PART_REAL
    !-----------------------------------------------
    !
    RO_S(1) = RHOS
    D_p0(1) = dia_phys
    
    DES_RADIUS(:) = ZERO
    PMASS(:) = ZERO
    PVOL(:) = ZERO
    OMOI(:) = ZERO
    RO_Sol(:) = ZERO 
    DES_POS_OLD(:,:) = ZERO
    DES_POS_NEW(:,:) = ZERO
    DES_VEL_OLD(:,:) = ZERO
    DES_VEL_NEW(:,:) = ZERO
    FC(:,:) = ZERO
    FN(:,:) = ZERO
    FT(:,:) = ZERO
    TOW(:,:) = ZERO
    OMEGA_OLD(:,:) = ZERO
    OMEGA_NEW(:,:) = ZERO

    FNS1(:) = ZERO
    FTS1(:) = ZERO

    NEIGHBOURS(:,:) = -1
    PN(:,:) = -1
    PV(:,:) = 1
    PFN(:,:,:) = ZERO
    PFT(:,:,:) = ZERO

    NEIGHBOURS(:,1) = 0
    PN(:,1) = 0
    PV(:,1) = 1
    PIJK(:,:) = ZERO

    DES_WALL_POS(:,:) = ZERO
    DES_WALL_VEL(:,:) = ZERO
    IS_MOBILE(1:PARTICLES) = .TRUE.
    CAUSE_MOTION(1:PARTICLES) = .TRUE.

    IF(GENER_CONFIG_CASE) THEN 
       if(I_AM_NODE_ZERO)WRITE(*,*) 'DES_INIT ARRAYS: GENER_CONFIG_CASE', GENER_CONFIG_CASE 
       DO LL = 1, PARTICLES 
          DES_POS_OLD(LL,1) = XC_GENER(LL,1)
          DES_POS_OLD(LL,2) = (XC_GENER(LL,2))
          DES_POS_OLD(LL,3) = (XC_GENER(LL,3))
          DES_RADIUS(LL) = RAD_GENER(LL)
          DES_VEL_OLD(LL,:) = ZERO
          OMEGA_OLD(LL,:) = ZERO
          DES_POS_NEW(LL,:) = DES_POS_OLD(LL,:)
          DES_VEL_NEW(LL,:) = DES_VEL_OLD(LL,:)
          OMEGA_NEW(LL,:) = OMEGA_OLD(LL,:)
       end DO
    ELSE
       
       if(I_AM_NODE_ZERO)WRITE(*,*) 'DES_INIT ARRAYS: GENER_CONFIG_CASE', GENER_CONFIG_CASE 
!!$       CALL init_random_seed
!!$       CALL random_number(TEST_PART_REAL)
!!$       Write(*,*)'TEST_PART_REAL = ', TEST_PART_REAL
!!$       TEST_PART = 1 + INT((PARTICLES-1)*TEST_PART_REAL)
!!$       TEST_PART = PARTICLES/2
		  TEST_PART = 1
       DO LL = 1, PARTICLES 
          DES_POS_OLD(LL,1) = (XC(LL,1)-one)*dx
          DES_POS_OLD(LL,2) = (XC(LL,2)-one)*dy
          DES_POS_OLD(LL,3) = (XC(LL,3)-one)*dz
          DES_RADIUS(LL) = RADBDY(LL)*dx
          DES_VEL_OLD(LL,:) = VELBDY(LL,:)
          OMEGA_OLD(LL,:) = ZERO
          DES_POS_NEW(LL,:) = DES_POS_OLD(LL,:)
          DES_VEL_NEW(LL,:) = DES_VEL_OLD(LL,:)
          OMEGA_NEW(LL,:) = OMEGA_OLD(LL,:)
       end DO

       if(CAGE_SIMUL.and.I_AM_NODE_ZERO)then
          WRITE(*,'(A)')'CAGE SIMULATION FLAG IS TURNED ON '
          WRITE(*,'(A30,2x,I6)')'INDEX OF THE TEST PARTICLE IS : ', TEST_PART
          WRITE(*,'(A30,3(2x,I6))')'CELL OF THE TEST PARTICLE IS : ', PIJK(TEST_PART,:)
          WRITE(*,'(A30,3(2x,g17.8))')'POSITION OF THE TEST PARTICLE IS : ', DES_POS_NEW(TEST_PART,:)
          WRITE(*,'(A30,3(2x,g17.8))')'POSITION OF THE TEST PARTICLE IS : ', DES_VEL_OLD(TEST_PART,:)
          IS_MOBILE(1:PARTICLES) = .FALSE.
          IS_MOBILE(TEST_PART) = .TRUE.
          CAUSE_MOTION(1:PARTICLES) = .FALSE.
          CAUSE_MOTION(TEST_PART) = .TRUE.
          !READ(*,*)
       end if
       if(CAGE_SIMUL)then
          DO LL = 1, PARTICLES
				 DES_VEL_NEW(LL,:) = DES_VEL_OLD(LL,:)
             if(LL.ne.TEST_PART)DES_VEL_OLD(LL,:) =  ZERO
          END DO
       end if
    end IF
    

    RETURN
  END SUBROUTINE DES_INIT_ARRAYS

  SUBROUTINE CFUPDATEOLD
    
    USE dem_mod
    
    IMPLICIT NONE
    
    INTEGER LL
    !     
    !---------------------------------------------------------------------
    
    DO LL = 1, PARTICLES
       DES_POS_OLD(LL,:) = DES_POS_NEW(LL,:)
       DES_VEL_OLD(LL,:) = DES_VEL_NEW(LL,:)
       OMEGA_OLD(LL,:) = OMEGA_NEW(LL,:)
    END DO
  END SUBROUTINE CFUPDATEOLD
  
  SUBROUTINE CFNEWVALUES(L)
    
    IMPLICIT NONE
    
    !DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
    
    INTEGER L, KK, K, NSPLIT, CHECK
    INTEGER IJK, I, J, KKK 
    DOUBLE PRECISION TEMPTIME, PASSTIME, D(DIMN), DIST,  V, rhat(dimn), rhat2, force_tmpy
    !     
      !---------------------------------------------------------------------
      !     Calculate new values
      !---------------------------------------------------------------------
      !     
      CHECK = 0 

      IF(shrink) THEN 
         rhat(2) = DES_POS_NEW(L,2) - YLENGTH/2.d0
         !IF(ABS(rhat(2)).GT.2.d0*MAX_RADIUS) then
         !PRINT*, 'SHRINKING'
         rhat2 = rhat(2)*rhat(2)
         rhat(2)= rhat(2)/ABS(rhat(2))
         
         FC(L, 2) = FC(L,2) - 980.d0*PMASS(L)*rhat(2)*(1.d0-exp(-4.d0*rhat2/(YLENGTH**2.d0)))
         !endif
      ENDIF
      

      IF(collision_type.eq."softsphere") THEN
         if(.not.gener_config_case)then
            if(L.eq.-1)WRITE(*,'(3(A,3(2x,g12.5),/))') 'FC = ', FC(L,:)/PMASS(L), 'FORCE = ',RHOF*FORCE(L,:)/PMASS(L),&
                 ' GRAV = ',(ONE-RHOF/RHOS)*GRAV(:)
            !test_force(:) = test_force(:) + FC(L,:)
            if(only_dem)then
               FC(L, :) = FC(L,:)/PMASS(L) !+ GRAV(:) 
            else
               test_force(:) = test_force(:) + FC(L,:)
               contact_force(L,:) =  FC(L,:)
               !FC(L,:) = FC(L,:)/PMASS(L)+ RHOF*FORCE(L,:)/PMASS(L) + (one-rhof/rhos)*GRAV(:) - frame_accln(:)
               FC(L,:) = FC(L,:)/PMASS(L)+ RHOF*(PRES(L,:) + VISC(L,:))/PMASS(L) - frame_accln(:) - mpg(:)*PVOL(L)/PMASS(L)

! COMMENTING THIS, BECAUSE GRAVITY IS CHANGED TO MPG. IF BOTH GRAVITY AND MPG ARE AVAILABLE, THEN UNCOMMENT THIS
!               FC(L,:) = FC(L,:) + (one-rhof/rhos)*GRAV(:)

            end if
         else
            FC(L, :) = FC(L,:)/PMASS(L) !+ GRAV(:) 
         end if
         
			if (IS_MOBILE(L)) then
		      DES_VEL_NEW(L,:) = FC(L,:)             
		      DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID
		      OMEGA_NEW(L,:) = OMEGA_OLD(L,:)! + TOW(L,:)*OMOI(L)*DTSOLID
			else
		      DES_VEL_NEW(L,:) = zero
		      OMEGA_NEW(L,:) = zero
			endif

         if(.not.gener_config_case)then
            if(L.eq.-1)then
               WRITE(*,'(4(A,3(2x,g12.5),/))')&
                    'FORCE = ', FORCE(L,:),&
                    'FORCE = ', FC(L,:),&
                    'OLD = ', DES_VEL_OLD(L,:),&
                    'NEW = ', DES_VEL_NEW(L,:)
               READ(*,*)
            end if
         end if
         FC(L,:) = ZERO
         TOW(L,:) = ZERO
         
      end IF
      
      if(IS_MOBILE(L))then
         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + DES_VEL_NEW(L,:)*DTSOLID 
      else
         DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
      end if
      
      IF(DES_POS_NEW(L,2).LT.ZERO.AND..NOT.DES_PERIODIC_WALLS) THEN
         if (I_AM_NODE_ZERO) PRINT*,'POSITION LE ZERO FOR L = ', L, DES_VEL_NEW(L,:), DES_POS_NEW(L,:)
      ENDIF
      !PRINT*,'grav = ', GRAV(2), pgrad(2)*PVOL(L)/PMASS(L), RO_sol(L)
!!$      IF(.NOT.DO_NSEARCH) THEN
!!$         D(:) = DES_POS_NEW(L,:) - PPOS(L,:)
!!$         DIST = SQRT(DES_DOTPRDCT(D,D))
!!$         
!!$         NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L)
!!$         IF(DIST.GE.NEIGHBOR_SEARCH_DIST) DO_NSEARCH = .TRUE.
!!$      END IF
      
      ! Chacking if the partcile moved more than a dia in a solid time step   
      D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
      DIST = SQRT(DES_DOTPRDCT(D,D))
      IF(DIST.GE.DES_RADIUS(L)) THEN
			if (I_AM_NODE_ZERO) then
		      PRINT *,'MOVEMENT UNDESIRED: PARTICLE', L
		      PRINT*, 'DES POS NEW : ', DES_POS_NEW(L,:)
		      PRINT*, 'DES POS OLD : ', DES_POS_OLD(L,:)
			endif
         STOP
      END IF
      
      ! Periodic treatment
      IF(DES_PERIODIC_WALLS) THEN
         IF(DES_PERIODIC_WALLS_X) THEN
            IF(DES_POS_NEW(L,1).GT.EX2) THEN
               DES_POS_NEW(L,1) = DES_POS_NEW(L,1) - (EX2 - WX1)
               PIJK(L,1) = 2
            ELSE IF(DES_POS_NEW(L,1).LT.WX1) THEN
               DES_POS_NEW(L,1) = DES_POS_NEW(L,1) + (EX2 - WX1)
               PIJK(L,1) = IMAX1
            END IF
         END IF
         IF(DES_PERIODIC_WALLS_Y) THEN
            IF(DES_POS_NEW(L,2).GT.TY2) THEN
               DES_POS_NEW(L,2) = DES_POS_NEW(L,2) - (TY2 - BY1)
               PIJK(L,2) = 2
            ELSE IF(DES_POS_NEW(L,2).LT.BY1) THEN
               DES_POS_NEW(L,2) = DES_POS_NEW(L,2) + (TY2 - BY1)
               PIJK(L,2) = JMAX1
            END IF
         END IF
         IF(DES_PERIODIC_WALLS_Z) THEN
            IF(DES_POS_NEW(L,3).GT.NZ2) THEN
               DES_POS_NEW(L,3) = DES_POS_NEW(L,3) - (NZ2 - SZ1)
               PIJK(L,3) = 2
            ELSE IF(DES_POS_NEW(L,3).LT.SZ1) THEN
               DES_POS_NEW(L,3) = DES_POS_NEW(L,3) + (NZ2 - SZ1)
               PIJK(L,3) = KMAX1
            END IF
         END IF
      END IF

      RETURN
    END SUBROUTINE CFNEWVALUES

    SUBROUTINE IBMUPDATE
      IMPLICIT NONE
      INTEGER M
      DO M=1,NBODY
         XC(M,1) = DES_POS_NEW(M,1)/dx + one
         XC(M,2) = DES_POS_NEW(M,2)/dy + one 
         XC(M,3) = DES_POS_NEW(M,3)/dz + one 
         VELBDY(M,:) = DES_VEL_NEW(M,:)
         ANGV(M,:,1) = OMEGA_NEW(M,:)
      END DO
    END SUBROUTINE IBMUPDATE
    
    SUBROUTINE  HARD_SPHERE_COLLISION
      IMPLICIT NONE 
      INTEGER :: LL,IP, idim, COLL_I(1), COLL_J,PNO
      REAL(prcn) :: tij, rij(dimn),vij(dimn),tempr, rmax(dimn),rijsq,rpijsq,bij, discr, vijsq,tij_min, tij_tmp
      !DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      REAL(prcn) :: COLLTIME(PARTICLES),tmin, DLENGTH(dimn)
      INTEGER :: PARTNER(PARTICLES),COUNT

      COLLTIME(:) = LARGE_NUMBER
      PARTNER(:) = PARTICLES
      DLENGTH(1) = XE(IMAX1) - XE(1)
      DLENGTH(2) = YN(JMAX1) - YN(1)
      DLENGTH(3) = ZT(KMAX1) - ZT(1)
      rmax(:) = DLENGTH(:)/two
      
      DO LL=1,PARTICLES
         IF(NEIGHBOURS(LL,1).GT.1)THEN
            DO IP = 2, NEIGHBOURS(LL,1)+1
               PNO = NEIGHBOURS(LL,IP)
               rijsq = zero ! separation
               vijsq = zero
!!$               PRINT*,'IN HS, PNO = ', PNO
               IF(PNO.GT.LL)THEN
                  do idim = 1, dimn
                     tempr = DES_POS_NEW(LL,idim) - DES_POS_NEW(PNO,idim) ! compute the separation
                     ! in each dimension
                     
                     if((ABS(tempr)).gt.rmax(idim)) then
                        if(tempr.lt.zero) then
                           tempr = tempr + DLENGTH(idim)
                        else
                           tempr = tempr - DLENGTH(idim)
                        end if
!!$                        PRINT*,'tempr =', tempr
                     end if
                     
                     rij(idim) = tempr
                     vij(idim) = DES_VEL_NEW(LL,idim)-DES_VEL_NEW(PNO,idim)
                     rijsq = rijsq + tempr**2.d0
                     vijsq = vijsq + vij(idim)**2.d0
                  end do
!!$                  PRINT*,'rij=',rij(1:dimn), 'vij=',vij(1:dimn)
!!$                  bij = DES_DOTPRDCT(rij,vij)
                  bij = rij(1)*vij(1)+rij(2)*vij(2)+rij(3)*vij(3)
!!$                  PRINT*,'bij1 = ', bij, rij, vij
!!$                  READ(*,*)
                  IF(bij.lt.zero) THEN
                     rpijsq = (DES_RADIUS(LL)+DES_RADIUS(PNO))**2.d0
                     discr = bij ** 2 - vijsq * ( rijsq - rpijsq )
                     
!!$                     PRINT*,'bijsq = ', bij**2.d0
!!$                     PRINT*,'vijsq = ', vijsq
!!$                     PRINT*,'rijsq = ', rijsq
!!$                     PRINT*,'rpijsq = ', rpijsq
!!$                     PRINT*,'discr =' ,discr
                     if ( discr .gt. 0.d0 ) then
                        
                        tij = ( -bij - sqrt ( discr ) ) / vijsq
                        if(tij.lt.colltime(LL))then
                           COLLTIME(LL) = tij
                           PARTNER(LL) = PNO
                        end if
                     end if
                  end IF
               END IF
            END DO
         END IF
      END DO
      
      
      tij_min = MINVAL(COLLTIME(1:PARTICLES))
      Write(*,'(A,2x,g17.6,2x,A,2x,g17.6)')'COLLISION TIME = ', tij_min, 'DTSOLID = ', dths
      COLL_I = MINLOC(COLLTIME(1:PARTICLES))
      COLL_J = PARTNER(COLL_I(1))
      Write(*,'(A,2(2x, I6))')'MINIMUM COLLISION TIME FOR PARTICLE Nos', COLL_I(1),COLL_J
      !PRINT*,'VELOCITY OF THE PARTICLES', DES_VEL_NEW(COLL_I(1),:),DES_VEL_NEW(COLL_J,:)
      !READ(*,*)
      ! If minimum collision time is less than the time remaining in hardsphere, then set DTSOLID to minimum collision time so that the particle positions are updated for this time

      IF(tij_min.lt.dths) DTSOLID = tij_min

      DO LL=1, PARTICLES
         CALL CFNEWVALUES(LL)
      END DO
      
      IF(tij_min.lt.dths)THEN
         COUNT = PARTICLES
         Do While(COUNT.gt.1)
            tmin = MINVAL(COLLTIME(1:COUNT))
            COLL_I = MINLOC(COLLTIME(1:COUNT))
            COLL_J = PARTNER(COLL_I(1))
!            PRINT*,'DISTANCE BETWEEN PARTICLES = ', DES_POS_NEW(COLL_I(1),1)-DES_POS_NEW(COLL_J,1)
!            READ(*,*)
            if(tmin.eq.tij_min)THEN
               PRINT*,'CALLING BUMP FOR PARTICLES', COLL_I(1),COLL_J
               !READ(*,*)
               CALL BUMP(COLL_I(1),COLL_J)
               tij_tmp = COLLTIME(COUNT)
               COLLTIME(count) = COLLTIME(coll_i(1))
               COLLTIME(coll_i(1)) = tij_tmp
               COUNT = COUNT-1
               IF(GENER_CONFIG_CASE)THEN
                  MINCOLLS(COLL_I(1)) = MINCOLLS(COLL_I(1))+1
                  MINCOLLS(COLL_J) = MINCOLLS(COLL_J)+1
               ELSE
                  TOTCOLLS = TOTCOLLS + 1
               END IF
            else
               goto 222
            end if
         end Do
      END IF
222   CONTINUE
      PRINT*,'DONE WITH HARD SPHERE COLLISION', dt, tij_min
    END SUBROUTINE HARD_SPHERE_COLLISION
    
    SUBROUTINE BUMP ( I, J )
      
      
         ! *******************************************************************
         ! ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
         ! **                                                               **
         ! ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
         ! *******************************************************************
      IMPLICIT NONE 

         INTEGER::     I, J, idim
         DOUBLE PRECISION :: m1, m2 

         double precision  rij(dimn), factor, tempr, rpijsq, normal(dimn)
         double precision  vij(dimn), delv(dimn)
         double precision  vijsq, rmax(dimn), DLENGTH(dimn)
         !DOUBLE PRECISION :: DES_DOTPRDCT 
         !     *******************************************************************
         m1 = PMASS(I)
         m2 = PMASS(J)
         vijsq = zero
         DLENGTH(1) = XE(IMAX1) - XE(1)
         DLENGTH(2) = YN(JMAX1) - YN(1)
         DLENGTH(3) = ZT(KMAX1) - ZT(1)
         rmax(:) = DLENGTH(:)/two
         if (I_AM_NODE_ZERO) PRINT*,'MASS = ', m1, m2
         do idim = 1, dimn
            tempr = DES_POS_NEW(I,idim) - DES_POS_NEW(J,idim) ! compute the separation
            ! in each dimension
            if((ABS(tempr)).gt.rmax(idim)) then
               if(tempr.lt.zero) then
                  tempr = tempr + DLENGTH(idim)
               else
                  tempr = tempr - DLENGTH(idim)
               end if
            end if
            rij(idim) = tempr
            vij(idim) = DES_VEL_NEW(I,idim)-DES_VEL_NEW(J,idim)

            vijsq = vijsq + vij(idim)**2.d0
         end do
         normal(:) = rij(:)/DSQRT(rij(1)**2.d0 + rij(2)**2.d0+rij(3)**2.d0)
!!$         PRINT*,'rij,vij = ', rij, vij
!!$         rpijsq = (DES_RADIUS(I)+DES_RADIUS(J))**2.d0
         factor = DOT_PRODUCT(normal(1:DIMN),vij(1:DIMN))
!!$         factor = factor/rpijsq
!!$         factor = DES_DOTPRDCT(rij,vij)/rpijsq
!!$         PRINT*, 'FACTOR=', factor
         IF(factor.GT.0.d0)  THEN 
            WRITE(*,*)'ACHTUNG'
            WRITE(*,*) 'factor in bump GT zero = ', factor
            !STOP
         ENDIF
         if (I_AM_NODE_ZERO) PRINT*,'coeff_rest = ', coeff_rest

         do idim=1,dimn
            delv(idim) = -factor * normal(idim)
            DES_VEL_NEW(i,idim) = DES_VEL_NEW(i,idim) + delv(idim) * m2 *(1.d0+ coeff_rest)/(m1+m2)
            DES_VEL_NEW(j,idim) = DES_VEL_NEW(j,idim) - delv(idim) * m1 *(1.d0+ coeff_rest)/(m1+m2)
         end do
         !PRINT*,'NEW VELOCITIES=',DES_VEL_NEW(i,:),DES_VEL_NEW(j,:)
         !READ(*,*)
         return
       END SUBROUTINE BUMP
            
  
       SUBROUTINE  SOFT_SPHERE_COLLISION

         IMPLICIT NONE 
         !-----------------------------------------------
         !     L o c a l   V a r i a b l e s
         !-----------------------------------------------
         !     
         INTEGER NN, LN, I, J, K, NSN, NP, L, IJK, GTC_COUNT, GTC_FACTOR, idim
         DOUBLE PRECISION TEMP_DTS, DTSOLIDTEMP , MEAN_VEL(ndim)
         !CHARACTER*5 FILENAME
         !     Logical to see whether this is the first entry to this routine
         LOGICAL,SAVE:: FIRST_PASS = .TRUE.
         LOGICAL DES_CONTINUUM_COUPLED_FT
         DOUBLE PRECISION  pgrad_tmp(1:DIMN), GRAV_TMP(1:DIMN), grantemp
         DOUBLE PRECISION:: xmintmp, xmax, ymin, ymax, Ax, bx, ay, by, az, bz , mean_free_path, t_mfp

         !     Force calculation         
         !IF(DES_PERIODIC_WALLS) THEN
         CALL CALC_FORCE_DES
         
         
         DO LN = 1, PARTICLES
            CALL CFNEWVALUES(LN)
         END DO
!!$         do ln = 1, ndim
!!$            WRITE(*,*) 'USMEAN = ', SUM(des_vel_new(1:nbody,ln))/real(nbody,prcn)
!!$         end do
         !CALL FIND_CELL_INDEX
         !CALL PARTICLES_IN_CELL
         
         !CALL GRID_BASED_NEIGHBOUR_SEARCH
         
         IF(GENER_CONFIG_CASE.AND.(MOD(IFAC,1000).EQ.0)) THEN

            do idim = 1, ndim
               mean_vel(idim) = SUM(DES_VEL_NEW(1:PARTICLES,idim))/real(particles,prcn)
            end do
            mean_vel = zero
            grantemp = zero
            do LN = 1, PARTICLES
               do idim = 1, ndim
                  grantemp = grantemp + half*PMASS(LN)*(DES_VEL_NEW(LN,idim)-mean_vel(idim))**2.d0
               end do
            end do
            grantemp = grantemp/real(particles,prcn)
            WRITE(*,*) 'MAX_OVERLAP, IFAC = ', OVERLAP_MAX, IFAC, grantemp
         end IF

         if(TRIM(input_type).eq.'lubtest')then
            WRITE(*,'(3(A,3(2x,g12.5),/))')&
                 'FORCE = ', FORCE(1,:),&
                 'OLD = ', DES_VEL_NEW(1,:),&
                 'NEW = ', DES_VEL_NEW(2,:)
            
            WRITE(*,*)'POSITIONS = ', DES_POS_NEW(2,:)-DES_POS_NEW(1,:)
            READ(*,*)
         end if


         IF(TEST_YMAXVAL) THEN

            IF((MOD(IFAC,1000).EQ.0)) THEN
               if (I_AM_NODE_ZERO) PRINT*,'MAX - MIN = ', MAXVAL(DES_POS_NEW(1:PARTICLES,2))-MINVAL(DES_POS_NEW(1:PARTICLES,2))

            end IF
            IF(IFAC.EQ.1) WRITE(*,*) 'TARGET MAX-MIN =', YMAXVAL-2.d0*MAX_RADIUS
            !WRITE(*,*)'DEKHO', MAXVAL(DES_POS_NEW(1:PARTICLES,2)),MINVAL(DES_POS_NEW(1:PARTICLES,2)),YMAXVAL, MAX_RADIUS*2.d0
            IF(MAXVAL(DES_POS_NEW(1:PARTICLES,2))-MINVAL(DES_POS_NEW(1:PARTICLES,2)).LT.(YMAXVAL-2.d0*MAX_RADIUS)) THEN

               WRITE(*,*) 'STOPPING THE SIMULATION BECAUSE MAX Y LT YMAXVAL'
               !WRITE(*,*) 'DES POS MAX  = ', MAXVAL(DES_POS_NEW(1:PARTICLES,2))
               !WRITE(*,*) 'DES POS MIN  = ', MINVAL(DES_POS_NEW(1:PARTICLES,2))
               WRITE(*,*) 'YMAXVAL  = ', YMAXVAL 
               ymin = MINVAL(DES_POS_NEW(1:PARTICLES,2)) - MAX_RADIUS
               ymax = MAXVAL(DES_POS_NEW(1:PARTICLES,2)) + MAX_RADIUS

               Ay = YMAXVAL/(ymax-ymin)
               By = -Ay*ymin

               DO L = 1, PARTICLES 

                  DES_POS_NEW(L,2) = Ay*DES_POS_NEW(L,2) + By
               end DO

               GOTO 100
            end IF
         end IF

         !     Write Restart for DEM only case

         S_TIME = S_TIME + DTSOLID

         !END DO                    ! end do NN = 1, FACTOR

         return 
100      continue
         SHRINK = .FALSE.
         if (I_AM_NODE_ZERO) PRINT*,'MAX - MIN = ', MAXVAL(DES_POS_NEW(1:PARTICLES,2))-MINVAL(DES_POS_NEW(1:PARTICLES,2))
         WRITE(*,*) 'TARGET MAX-MIN =', YMAXVAL-2.d0*MAX_RADIUS
         RETURN
         !IF(TEST_YMAXVAL) THEN
         !ENDIF

       END SUBROUTINE SOFT_SPHERE_COLLISION


    SUBROUTINE init_particles_jn
      USE randomno
      USE dem_mod
      IMPLICIT NONE 
      INTEGER :: i,j, k , ip, idim
      REAL*8 :: umf0(dimn), rsf(DIMN, DIMN), rhat(3), mean_vel(ndim), grantemp

      WRITE(*,*) 'INITIALIZING NORMAL VELOCITY DISTRIBUTION'
      WRITE(*,*) 'MEAN  = ', ZERO, ' AND VARIANCE = ', pvel_var
      do j=1,DIMN
         umf0(j)=zero
         do i=1,dimn
            if(i.eq.j)then
               rsf(i,j)=pvel_var/three
            else
               rsf(i,j)=0.0
            endif
         enddo
      enddo
      CALL jn_dist(DES_VEL_OLD(1:PARTICLES, 1:DIMN),particles,dimn,umf0,rsf)
      DO ip = 1, particles 

         IF(SHRINK) THEN 
            rhat(2) = DES_POS_NEW(ip,2) - YLENGTH/2.d0

            rhat(2)= rhat(2)/ABS(rhat(2))

            DES_VEL_OLD(ip,2) = ABS(DES_VEL_OLD(ip,2))*rhat(2) !Direct the y- velocity inwards 
         ENDIF

         DES_VEL_NEW(ip,:) = DES_VEL_OLD(ip,:)
      end DO
      grantemp = zero
      mean_vel = zero
      do IP = 1, PARTICLES
         do idim = 1, ndim
            grantemp = grantemp + half*PMASS(IP)*(DES_VEL_NEW(IP,idim)-mean_vel(idim))**2.d0
         end do
      end do
      grantemp = grantemp/real(particles,prcn)
      !WRITE(*,*) 'TOTAL INITIAL KE = ', grantemp
      !READ(*,*)
    end SUBROUTINE init_particles_jn

    SUBROUTINE init_random_seed
      Implicit None

      INTEGER              :: isize,idate(8)
      INTEGER,ALLOCATABLE  :: iseed(:)
      
      CALL DATE_AND_TIME(VALUES=idate)
      CALL RANDOM_SEED(SIZE=isize)
      ALLOCATE( iseed(isize) )
      CALL RANDOM_SEED(GET=iseed)
      iseed = iseed * (idate(8)-500) ! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)
      
      DEALLOCATE( iseed )
      
    END SUBROUTINE init_random_seed


  end MODULE collision_mod
  
