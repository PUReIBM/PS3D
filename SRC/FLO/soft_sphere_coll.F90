!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Purpose: DES calculations of force acting on a particle,            C
!           its velocity and its position                              C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: Now includes particle-wall interaction history.           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
MODULE soft_spring
  USE precision
  USE constants
  USE global_data
  USE dem_mod
  !PRIVATE
  LOGICAL :: PARTICLE_SLIDE
  
CONTAINS 
  SUBROUTINE CALC_FORCE_DES
    
    
    IMPLICIT NONE
    
    INTEGER LL, I, II, K, LC, IW, KK, TEMP, TEMPN, JJ, J, FOCUS_PART2 , NEIGH_L, NEIGH_I
    INTEGER NI, IJ, WALLCHECK, NLIM, N_NOCON, IP2, FOCUS_PARTICLE, NI_I
    INTEGER KM1, KP1, IM1, IP1, JM1, JP1, PNO, NPG, PC(3),  OVERLAP_MAXP
    DOUBLE PRECISION OVERLAP_N, OVERLAP_T, TEMPX, TEMPY, TEMPZ, TEMPD
    DOUBLE PRECISION V_REL_TRANS(DIMN), V_SLIP(DIMN)
    DOUBLE PRECISION NORMAL(DIMN), TANGENT(DIMN), DISTMOD
    DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, VSLIPMAG, MASS_I, MASS_J, MASS_EFF
    DOUBLE PRECISION TEMPFN(DIMN), TEMPFT(DIMN), DIST(DIMN), R_LM
    LOGICAL ALREADY_EXISTS
    LOGICAL CHECK_CON, ALREADY_NEIGHBOURS, OVERLAP_MAX_WALL 
    !     
                  !---------------------------------------------------------------------
!     Calculate new values
!---------------------------------------------------------------------
      OVERLAP_MAXP = UNDEFINED_I
      OVERLAP_MAX_WALL = .FALSE.
                          !     
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1
      FOCUS_PARTICLE = 0
      FOCUS_PART2 = 200000000
      IF (S_TIME.LE.DTSOLID) THEN
         V_REL_TRANS(:) = ZERO
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FC(:,:) = ZERO
         FN(:,:) = ZERO
         FT(:,:) = ZERO
!        PN(:,1) = 0
!        PN(:,2:MAXNEIGHBORS) = -1
!        PV(:,:) = 1
!         PFN(:,:,:) = ZERO
!         PFT(:,:,:) = ZERO
      END IF


      !DO LL = 1, PARTICLES
      !   IF(ABS(FC(LL,2)).GT.ZERO) THEN 
      !      PRINT*,'FORCE 2 GT 0', ABS(FC(LL,2))/PMASS(LL)
      !   end IF
      !ENDDo

!     
!---------------------------------------------------------------------
!     Calculate contact force and torque
!---------------------------------------------------------------------
!     
      DO LL = 1, PARTICLES
         
!         IF(LL.EQ.FOCUS_PARTICLE) THEN 
!            PRINT*, DES_POS_NEW(LL,1), DES_POS_NEW(LL,2)
!            PRINT*, DES_VEL_NEW(LL,1), DES_VEL_NEW(LL,2)
            !PRINT*, 'NEIGHS  = ',  PN(LL,:)
!         end IF
         
         NI = 0
         TEMPFN(:) = ZERO
         TEMPFT(:) = ZERO
         K = 0
         KK = 0
         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
            DO NI = 2, PN(LL,1)+1
               IF(PV(LL,NI-N_NOCON).EQ.0.AND.NI-N_NOCON.GE.2) THEN
                  PN(LL,NI-N_NOCON:NLIM-1) = PN(LL,NI-N_NOCON+1:NLIM) 
                  PV(LL,NI-N_NOCON:NLIM-1) = PV(LL,NI-N_NOCON+1:NLIM) 
                  PFN(LL,NI-N_NOCON:NLIM-1,:) = PFN(LL,NI-N_NOCON+1:NLIM,:) 
                  PFT(LL,NI-N_NOCON:NLIM-1,:) = PFT(LL,NI-N_NOCON+1:NLIM,:) 
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
                  NLIM = PN(LL,1) + 1
               ENDIF
            END DO
         ENDIF
         
!     Initializing rest of the neighbor list which is not in contact

         NLIM = MAX(2,PN(LL,1) + 2) 
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFN(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         
         NEIGH_MAX = MAX(NEIGH_MAX, PN(LL,1)+1)
         !     Initializing the neighbor list contact information when particles are not in contact

         IF (PN(LL,1).EQ.0) THEN
           PFN(LL,:,:) = ZERO
           PFT(LL,:,:) = ZERO
         END IF
         
!     Initializing the particle
         DO K = 2, MAXNEIGHBORS
            PV(LL,K) = 0
         END DO
         
         IF(WALLDTSPLIT) THEN
            WALLCHECK = 0
            DO IW = 1, NWALLS
               WALLCONTACT = 0
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
               IF(WALLCONTACT.EQ.1) THEN
                  WALLCHECK = 1
                  !PRINT*, 'WAll = ', IW, LL
                  I = PARTICLES + IW
                  
                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN
                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        end IF
                        
                     ENDDO
                  end IF
                  
                  CALL CFWALLPOSVEL(LL, IW)
                  DES_POS_NEW(I,:) = DES_WALL_POS(IW,:)
                  DES_VEL_NEW(I,:) = DES_WALL_VEL(IW,:)
                  OMEGA_NEW(I,:) = ZERO
                  DES_RADIUS(I) = DES_RADIUS(LL)
                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))

                  IF(R_LM - DISTMOD.gt.SMALL_NUMBER) then 
                     
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THen
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                        OVERLAP_MAX_WALL = .TRUE.
                     ENDIF
                     CHECk_CON = .TRUE.
                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        PRINT *,'DISTMOD IS ZERO', I,LL
                        STOP
                     END IF
                     CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                     
                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_N =  R_LM-DISTMOD!V_REL_TRANS_NORM*DTSOLID
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        OVERLAP_N = R_LM-DISTMOD!V_REL_TRANS_NORM*DTSOLID! 
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ENDIF
                  ELSE
                     CHECk_CON = .FALSE.
                     GOTO 200
                  endif
                  
                              
                  MASS_I = (PI*((2.d0*DES_RADIUS(LL))**3.d0)*RO_Sol(LL))/6.d0
                  MASS_J = MASS_I
                  MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)
                  ETA_N_W = 2.D0*SQRT(KN_W*MASS_EFF)*ABS(LOG(REAL_EN_WALL(1)))
                  
                  ETA_N_W = ETA_N_W/SQRT(PI*PI + (LOG(REAL_EN_WALL(1)))**2.0)
                  !DES_ETAT(I,J) = HALF*DES_ETAN(I,J)
                  ETA_T_W = HALF*ETA_N_W
                  
                  
                  
                  
                  FNS1(:) = -KN_W*((OVERLAP_N))*NORMAL(:)
                  FNS2(:) = -ETA_N_W*V_REL_TRANS_NORM*NORMAL(:)
                  
                  FTS1(:) = -KT_W*((OVERLAP_T)) *TANGENT(:)
                  FTS2(:) = -ETA_T_W*V_REL_TRANS_TANG*TANGENT(:)
                  
                  
                  FT(LL,:) = FTS1(:) + FTS2(:) 
                  FN(LL,:) = FNS1(:) + FNS2(:) 
                  
                  FN(LL, :) = FN(LL,:)! + PFN(LL,NI,:)
                  TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                  
                  CALL CFSLIDEWALL(LL, TANGENT, TEMPFT)
                  
                  CALL CFFCTOWALL(LL, NORMAL)

                  PFN(LL,NI,:) =  PFN(LL,NI,:) + FNS1(:)
                  
                  IF(.NOT.PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                  ELSE
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                     PARTICLE_SLIDE = .FALSE.
                  END IF
                  
                  IF(LL.EQ.FOCUS_PARTICLE.and.(IW.Eq.4.OR.IW.EQ.3)) THEN 
                     PRINT*,'                             '
                     PRINT*,'WAL CONTACT ON WITH WALL', IW
                     PRINT*, 'ALREADY_NEIGHBOURS? = ', ALREADY_NEIGHBOURS  
                     PRINT*, 'PV = ',PV(LL,:)
                     PRINT*,'PVEL = ;', DES_VEL_NEW(LL,:)
                     PRINT*,'OVERLAPS = ;', OVERLAP_N, OVERLAP_T
                     !PRINT*,' ETA_N = ', KN_W, ETA_N_W, KT_W, ETA_T_W
                     PRINT*, 'NORMAL = ', NORMAL(:), KN_W, ETA_N_W
                     PRINT*,'FN = ',FN(LL,:)
                     PRINT*,'FNS1 and FNS2 = ',FNS1(:), FNS2(:)
                     PRINT*,'PFN = ',PFN(LL,NI,:)
                     PRINT*,'TANGENT = ', TANGENT(:)
                     PRINT*,'PFT = ', PFT(LL,NI,:)
                     PRINT*,'                             '
                     !READ(*,*)
                  ENDIF
                  

               end IF!Wall Contact
200            continue
            end DO
         end IF !if(walldtsplit)
      
         
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)

               IF(I.GT.LL) THEN
                  
                  ALREADY_NEIGHBOURS=.FALSE.
                     
                  IF(PN(LL,1).GT.0) THEN
                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        end IF
                     ENDDO
                     
                     IF(ALREADY_NEIGHBOURS) THEN 
                        DO NEIGH_I = 2, PN(I,1)+1
                           IF(LL.EQ. PN(I,NEIGH_I)) THEN 
                              NI_I = NEIGH_I
                              EXIT
                           end IF
                        ENDDO
                     end IF
                     
                  end IF
                  
                  IF(DES_PERIODIC_WALLS) THEN
                     TEMPX = DES_POS_NEW(I,1)
                     TEMPY = DES_POS_NEW(I,2)
                     IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(I,3)        
                     TEMPD = ABS(DES_POS_NEW(LL,1) - DES_POS_NEW(I,1))
                     IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_X) THEN 
                        !PRINT*,'II, L = ', I,LL
                        !PRINT*,'OLD POS, I', DES_POS_NEW(I,:)
                        !PRINT*,'OLD POS, LL = ', DES_POS_NEW(LL,:)
                        IF(TEMPX.GT.DES_POS_NEW(LL,1)) THEN 
                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) - (EX2-WX1)
                         !  PRINT*,'NEW POS WEST= ', DES_POS_NEW(I,1)
                        ELSE
                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) + EX2 - WX1
                           
                         !  PRINT*,'NEW POS EAST = ', DES_POS_NEW(I,1)
                        ENDIF
                     ENDIF
         
                     TEMPD = ABS(DES_POS_NEW(LL,2) - DES_POS_NEW(I,2))
                     IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Y) THEN
                        IF(TEMPY.GT.DES_POS_NEW(LL,2)) THEN 
                           DES_POS_NEW(I,2) = DES_POS_NEW(I,2) - (TY2-BY1)
                        ELSE
                           DES_POS_NEW(I,2) = DES_POS_NEW(I,2) + (TY2-BY1)
                        ENDIF
                     ENDIF
                     
                     IF(DIMN.EQ.3) THEN
                        TEMPD = ABS(DES_POS_NEW(LL,3) - DES_POS_NEW(I,3))
                        IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Z) THEN 
                           IF(TEMPZ.GT.DES_POS_NEW(LL,3)) THEN 
                              DES_POS_NEW(I,3) = DES_POS_NEW(I,3) -(NZ2 - SZ1)
                           ELSE
                              DES_POS_NEW(I,3) = DES_POS_NEW(I,3) + (NZ2-SZ1)
                           ENDIF
                        ENDIF
                     ENDIF
                  END IF
                  
                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
                  
                  IF(DES_PERIODIC_WALLS) THEN
                     DES_POS_NEW(I,1) = TEMPX
                     DES_POS_NEW(I,2) = TEMPY
                     IF (DIMN.EQ.3) DES_POS_NEW(I,3) = TEMPZ              
                  END IF

                  IF(R_LM - DISTMOD.gt.SMALL_NUMBER) then 
                     IF(LL.EQ.FOCUS_PARTICLE) Print*, 'NEIGHBORS', NEIGHBOURS(LL,:)

                     
                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THen
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                        OVERLAP_MAX_WALL = .FALSE.
                     ENDIF

                     CHECk_CON = .TRUE.
                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        PRINT *,'DISTMOD IS ZERO FOR PART. PAIR', I,LL
                        STOP
                     END IF
                     CALL CFRELVEL(LL, I, V_REL_TRANS, V_REL_TRANS_NORM, V_REL_TRANS_TANG, TANGENT, NORMAL)
                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        PV(I,NI_I) = 1
                        OVERLAP_N =  R_LM-DISTMOD!V_REL_TRANS_NORM*DTSOLID
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        !pRINT*,'INITIATING CONTACT'
                        !PRINT*,'OVERLAP = ', ((R_LM-DISTMOD)/R_LM)*100.d0
                                !READ(*,*)
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        
                        PN(I,1) = PN(I,1) + 1
                        NI = PN(I,1) + 1
                        PN(I,NI) = LL
                        PV(I,NI) = 1

                        OVERLAP_N = R_LM-DISTMOD!V_REL_TRANS_NORM*DTSOLID! 
                        OVERLAP_T = ZERO!VRT*DTSOLID
                        
                     END IF
                  ELSE
                     GOTO 300
                  END IF
                            
                  MASS_I = (PI*((2.d0*DES_RADIUS(LL))**3.d0)*RO_Sol(LL))/6.d0
                  MASS_J = (PI*((2.d0*DES_RADIUS(I))**3.d0)*RO_Sol(I))/6.d0
                  MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)
                  ETA_DES_N = 2.D0*SQRT(KN*MASS_EFF)*ABS(LOG(REAL_EN(1,1)))

                  ETA_DES_N = ETA_DES_N/SQRT(PI*PI + (LOG(REAL_EN(1,1)))**2.0)
                  !DES_ETAT(I,J) = HALF*DES_ETAN(I,J)
                  ETA_DES_T = HALF*ETA_DES_N

                  FNS1(:) = -KN*((OVERLAP_N))*NORMAL(:)
                  FNS2(:) = -ETA_DES_N*V_REL_TRANS_NORM*NORMAL(:)
                           
                  FTS1(:) = -KT*((OVERLAP_T)) *TANGENT(:)
                  FTS2(:) = -ETA_DES_T*V_REL_TRANS_TANG*TANGENT(:)
                  
                  
                  FT(LL,:) = FTS1(:) + FTS2(:) 
                  FN(LL,:) = FNS1(:) + FNS2(:) 
                  
                  FN(LL, :) = FN(LL,:)! + PFN(LL,NI,:)
                  TEMPFT(:) = FT(LL, :) + PFT(LL,NI,:)
                  
                  IF(LL.EQ.FOCUS_PARTICLE) THEN 
                     PRINT*, 'I = ', I
                     PRINT*,'ETAs =  ', ETA_DES_N, ETA_DES_T
                     PRINT*,'overlap = ', overlap_n, (R_LM - DISTMOD)*100.d0/R_LM
                                      
                     PRINT*,'rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, 'FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, 'PFN = ', PFN(LL,NI,:)
                     PRINT*, 'PFT = ', PFT(LL,NI,:)
                     PRINT*, 'FORCEST = ', FT(LL,:)
                     PRINT*, 'FORCESN = ', FN(LL,:)
                     PRINT*, 'FORCEST = ', FT(LL,:)
                     READ(*,*)
                  ENDIF
                  CALL CFSLIDE(LL, TANGENT, TEMPFT)
                  
                  PFN(LL,NI,:) = PFN(LL,NI,:) +  FNS1(:)

                  CALL CFFCTOW(LL, I, NORMAL)

                  IF(.NOT.PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                  ELSE
                     PFT(LL,NI,:) = PFT(LL,NI,:) + FTS1(:)
                     PARTICLE_SLIDE = .FALSE.
                  END IF
                  
                  !!impulse is effectively doubled for wall interactions
                  IF(LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1)&
                        ,'FNy=',FN(LL,2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1)&
                       ,'FNy=',FN(LL,2)
                     END IF
                     CLOSE (1)
                  END IF
                  !--   ND DEBUGGING
                  
                  IF(LL.EQ.FOCUS_PARTICLE) Print*, 'PN', PN(LL,:)
               END IF

               
300            CONTINUE
               
               
            END DO              !II = 2, NEIGHBOURS(LL,1)+I
         END IF                 !(NEIGHBOURS(LL,1).GT.0)
         
         IF((NEIGHBOURS(LL,1).EQ.0).AND.(WALLCHECK.EQ.0)) THEN
            FN(LL,:) = ZERO
            FNS1(:) = ZERO
            FT(LL,:) = ZERO
            FTS1(:) = ZERO
            FC(LL,:) = ZERO
            TOW(LL,:) = ZERO
            OMEGA_NEW(LL,:) = ZERO
            
         END IF

      END DO
      
      IF(CAGE_SIMUL.and.(.not.GENER_CONFIG_CASE))then
         DO LL = 1, PARTICLES
!!$            if(.not.(IS_MOBILE(LL)))then
!!$               DO NEIGH_L = 2, PN(LL,1)+1  ! Loop over contact neighbours
!!$                  NI = PN(LL,NEIGH_L)
!!$                  if(CAUSE_MOTION(NI))then
!!$                     IS_MOBILE(LL) = .TRUE. ! If any of the contact neighbourse can cause motion, move this particle
!!$                     Write(*,*)' No. of CONTACTS FOR PARTICLE ', LL, ' = ', PN(LL,1)
!!$                     Write(*,*)'PARTICLE ', LL, ' IN CONTACT WITH DRIVER PARTICLE ', NI, 'MOVING THIS PARTICLE.'
!!$                     Write(*,*)'INDEX OF TEST PARTICLE ', TEST_PART
!!$                     EXIT
!!$                  end if
!!$               END DO
!!$            end if
!!$
            !if(LL.eq.test_part)then
               !write(*,*)'PN = ', PN(LL,1), PN(LL,2:PN(LL,1)+1)
               !READ(*,*)
            !end if
            if(LL.ne.TEST_PART)then
               IS_MOBILE(LL) = .FALSE. 
               
               DO NEIGH_L = 2, PN(LL,1)+1 ! Loop over contact neighbours
                  
                  NI = PN(LL,NEIGH_L)
                  !Write(*,*)'I AM HERE :', PN(LL,1)
                  !READ(*,*)
                  if(CAUSE_MOTION(NI))then
                     IS_MOBILE(LL) = .TRUE. ! If any of the contact neighbourse can cause motion, move this particle
                     Write(*,*)' No. of CONTACTS FOR PARTICLE ', LL, ' = ', PN(LL,1)
                     Write(*,*)'PARTICLE ', LL, ' IN CONTACT WITH DRIVER PARTICLE ', NI, 'MOVING THIS PARTICLE.'
                     Write(*,*)'INDEX OF TEST PARTICLE ', TEST_PART
                     EXIT
                  end if
               END DO
            end if
         END DO
      END IF
      !WRITE(*,*) 'MAX NEIGHBORS  = ', NEIGH_MAX
      !WRITE(*,*) 'MAX OVERLAP %  = ', OVERLAP_MAX

!-------------------------------------------------------------------
!     Update old values with new values
!-------------------------------------------------------------------
      !WRITE(*,*) 'PART MAX OVERLAPPED  AND WITH WALL ? = ',OVERLAP_MAXP, OVERLAP_MAX_WALL
      !   WRITE(*,*) DES_POS_NEW(1), DES_POS_NEW(2), ABS(DES_POS_NEW(1)- (DES_POS_NEW(2))/(DES_RADIUS(1)+DES_RADIUS(2)))*100
      !   WRITE(*,*) 'neigh 1 = ', NEIGHBOURS(1, 2:NEIGHBOURS(1,1)+1)
      !   IF(OVERLAP_MAX.GT.ZERO) THEN
      !      READ(*,*)
            
      !   ENDIF
      RETURN
      END SUBROUTINE CALC_FORCE_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLCONTACT(WALL, L, WALLCONTACTI)            C
!  Purpose: DES - Checking for contact with walls                      C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLCONTACT(WALL, L, WALLCONTACTI)

      Use dem_mod
      IMPLICIT NONE
      
      INTEGER L, I, K, WALL, WALLCONTACTI	
      DOUBLE PRECISION A, OMEGA, OOMEGA2, ASINOMEGAT 
      
!     
!---------------------------------------------------------------------
!     Checking if a particle is in contact with any of the walls
!---------------------------------------------------------------------
!     


      A = ZERO
      OMEGA = ZERO
      ASINOMEGAT = ZERO
      IF(DES_F.NE.ZERO) THEN
         OMEGA = 2.0D0*PI*DES_F
         OOMEGA2 = ONE/(OMEGA**2)
         A = DES_GAMMA*GRAV(2)*OOMEGA2
         ASINOMEGAT = A*SIN(OMEGA*S_TIME)
      END IF

      WALLCONTACTI = 0

      IF(WALL.EQ.1.AND.(.NOT.DES_PERIODIC_WALLS_X)) THEN
         IF((DES_POS_NEW(L,1)-WX1).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.2.AND.(.NOT.DES_PERIODIC_WALLS_X)) THEN
         IF((EX2-DES_POS_NEW(L,1)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.3.AND.(.NOT.DES_PERIODIC_WALLS_Y)) THEN
         IF((DES_POS_NEW(L,2)-(BY1+ASINOMEGAT)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.4.AND.(.NOT.DES_PERIODIC_WALLS_Y)) THEN
         IF((TY2-DES_POS_NEW(L,2)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.5.AND.(.NOT.DES_PERIODIC_WALLS_Z)) THEN
         IF((DES_POS_NEW(L,3)-SZ1).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.6.AND.(.NOT.DES_PERIODIC_WALLS_Z)) THEN
         IF((NZ2-DES_POS_NEW(L,3)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF
      END IF

      RETURN
    END SUBROUTINE CFWALLCONTACT
    


!                                                                      C
!  Module name: CFWALLPOSVEL(L, I)
!  Purpose:  DES -Calculate the position and velocity of wall particle C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLPOSVEL(L, I)

      IMPLICIT NONE

      INTEGER L, I 
      DOUBLE PRECISION A, OMEGA_W, F, COSOMEGAT, SINOMEGAT
      DOUBLE PRECISION DES_R, OOMEGAW2
      
!     
!---------------------------------------------------------------------
!     Assigning wall position and velocity
!---------------------------------------------------------------------
!     

      A = ZERO
      OMEGA_W = ZERO
      IF(DES_F.NE.ZERO) THEN
         OMEGA_W = 2.0d0*Pi*DES_F
         OOMEGAW2 = ONE/(OMEGA_W**2)
         A = DES_GAMMA*GRAV(2)*OOMEGAW2
         SINOMEGAT = SIN(OMEGA_W*S_TIME)
         COSOMEGAT = COS(OMEGA_W*S_TIME)
      END IF

      DES_R = DES_RADIUS(L)
      IF(WALLREFLECT) THEN
         DES_R = ZERO
      END IF

      DES_WALL_VEL(I,1) = ZERO
      DES_WALL_VEL(I,2) = ZERO
      IF(DIMN.EQ.3) DES_WALL_VEL(I,3) = ZERO

      IF(I.EQ.1) THEN
         DES_WALL_POS(I,1) = WX1 - DES_R	
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         WALL_NORMAL(1,1) = -ONE
         WALL_NORMAL(1,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(1,3) = ZERO

      ELSE IF(I.EQ.2) THEN
         DES_WALL_POS(I,1) = EX2 + DES_R
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         WALL_NORMAL(2,1) = ONE
         WALL_NORMAL(2,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(2,3) = ZERO

      ELSE IF(I.EQ.3) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = BY1 - DES_R + (A*SINOMEGAT)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         DES_WALL_VEL(I,1) = lid_vel
         DES_WALL_VEL(I,2) =  A*OMEGA_W*COSOMEGAT
         WALL_NORMAL(3,1) = ZERO
         WALL_NORMAL(3,2) = -ONE
         IF(DIMN.EQ.3) WALL_NORMAL(3,3) = ZERO

      ELSE IF(I.EQ.4) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = TY2 + DES_R
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         DES_WALL_VEL(I,1) = lid_vel
         WALL_NORMAL(4,1) = ZERO
         WALL_NORMAL(4,2) = ONE
         IF(DIMN.EQ.3) WALL_NORMAL(4,3) = ZERO

      ELSE IF(I.EQ.5) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = SZ1 - DES_R
         WALL_NORMAL(5,1) = ZERO
         WALL_NORMAL(5,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(5,3) = -ONE

      ELSE IF(I.EQ.6) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = NZ2 + DES_R
         WALL_NORMAL(6,1) = ZERO
         WALL_NORMAL(6,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(6,3) = ONE
      END IF
      
      RETURN
      END SUBROUTINE CFWALLPOSVEL



      SUBROUTINE CFRELVEL(L, II, VRELTRANS,VRN, VRT,  TANGNT, NORM)
      
      IMPLICIT NONE     
      !DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      
      INTEGER L, KK, II
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION TANMOD, TANGNT(DIMN), NORM(DIMN), VSLIP(DIMN), V_ROT(DIMN), OMEGA_SUM(DIMN), VRN, VRT
!
!-----------------------------------------------------------------------

      
      V_ROT(:) = ZERO
      
      
      VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DES_RADIUS(L)+ OMEGA_NEW(II,:)*DES_RADIUS(II)
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DES_RADIUS(L)+ OMEGA_NEW(II,1)*DES_RADIUS(II)
         OMEGA_SUM(2) = ZERO
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)
      
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

      VRN = DES_DOTPRDCT(VRELTRANS,NORM)
      
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
      TANMOD = ZERO
      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.0.) THEN
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      END IF
      
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)

      RETURN
      END SUBROUTINE CFRELVEL
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIDEWALL(L, TANGNT)                                 C
!  Purpose: DES - Calculate slide between particles and walls          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIDEWALL(L, TANGNT, TEMP_FT)
      
      IMPLICIT NONE

      !DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TEMP_FT(DIMN), TEMP_FN(DIMN)
!     
!---------------------------------------------------------------------


      !TEMP_FT(:) = FTS1(:)
      TEMP_FN(:) = FN(L, :)

      FTMD = SQRT(DES_DOTPRDCT(TEMP_FT,TEMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TEMP_FN,TEMP_FN))

      IF (FTMD.GT.(MEW*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         FT(L,:) = - MEW_W*FNMD*TANGNT(:)
      ELSE
         FT(L,:) = TEMP_FT(:)
      END IF
      

      RETURN
      END SUBROUTINE CFSLIDEWALL



      SUBROUTINE CFFCTOWALL(L, NORM)
     
      IMPLICIT NONE
      
      INTEGER L, K
      DOUBLE PRECISION NORM(DIMN), CROSSP(DIMN), FT1(DIMN)
!---------------------------------------------------------------------
!     

        FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 

        FT1(:) = FT(L,:)

        
         IF(DIMN.EQ.3) THEN 
            CALL DES_CROSSPRDCT(CROSSP, NORM, FT1)
            TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*CROSSP(:)
         ELSE 
            CROSSP(1) = NORM(1)*FT1(2) - NORM(2)*FT1(1)
            TOW(L,1) = TOW(L,1) + DES_RADIUS(L)*CROSSP(1)
         endif 
         !CALL DES_CROSSPRDCT(CROSSP, NORM, FT1)         

         !TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*CROSSP(:)
         
         !IF(L.EQ.29) PRINT*,'FOR = ', L, FC(L,:)
      RETURN
      END SUBROUTINE CFFCTOWALL

      SUBROUTINE CFSLIDE(L, TANGNT, TEMP_FT)

      IMPLICIT NONE

      !DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
      INTEGER L, K
      DOUBLE PRECISION FTMD, FNMD, TANGNT(DIMN)
      DOUBLE PRECISION TEMP_FT(DIMN), TEMP_FN(DIMN)
!     
!---------------------------------------------------------------------

      TEMP_FN(:) = FN(L, :)

      FTMD = SQRT(DES_DOTPRDCT(TEMP_FT,TEMP_FT))
      FNMD = SQRT(DES_DOTPRDCT(TEMP_FN,TEMP_FN))

      IF (FTMD.GT.(MEW*FNMD)) THEN
          PARTICLE_SLIDE = .TRUE.
          FT(L,:) = - MEW*FNMD*TANGNT(:)
      ELSE
         FT(L, :) = TEMP_FT(:)
      END IF

      RETURN
      END SUBROUTINE CFSLIDE


      SUBROUTINE CFFCTOW(L, II,  NORM)
     
      IMPLICIT NONE
      
      INTEGER L, II, K
      DOUBLE PRECISION NORM(DIMN), CROSSP(DIMN), FT1(DIMN), FT2(DIMN)
!---------------------------------------------------------------------
!     


         FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 
         FC(II,:) = FC(II,:) - FN(L,:) - FT(L,:)

         FT1(:) = FT(L,:)

         IF(DIMN.EQ.3) THEN 
            CALL DES_CROSSPRDCT(CROSSP, NORM, FT1)
            TOW(L,:) = TOW(L,:) + DES_RADIUS(L)*CROSSP(:)
            TOW(II,:) = TOW(II,:) + DES_RADIUS(II)*CROSSP(:)
	!Remember the torque is r cross F_T, which, compared to i particle, are both negative for the j particl.
	!Therefore, the toqrue, unlike tangential and normal contact forces, will be in the same direction for 
	!both the particles making the pair 
         ELSE 
            CROSSP(1) = NORM(1)*FT1(2) - NORM(2)*FT1(1)
            TOW(L,1) = TOW(L,1) + DES_RADIUS(L)*CROSSP(1)
            TOW(II,1) = TOW(II,1) + DES_RADIUS(II)*CROSSP(1)
         endif 

	!         IF(L.EQ.29) PRINT*,'FOR = ', L, II, FC(L,:)

      RETURN
      END SUBROUTINE CFFCTOW




    end MODULE SOFT_SPRING
    
