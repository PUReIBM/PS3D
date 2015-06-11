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

MODULE post_static
#include "../FLO/ibm.h"
  !
  ! AUTHOR: Sudheer Tenneti
  ! POST PROCESSOR CODE FOR MEAN DRAG MODELS

  USE precision  
  USE restart_funcs
  USE general_funcs
  USE postproc_funcs
  USE constants  
  USE randomno
  USE useful_functions
  USE dependent_functions
  USE nlcalc
  USE post_global_data

  IMPLICIT NONE
 
  PRIVATE
  REAL(prcn), ALLOCATABLE, DIMENSION(:) :: RAD_PHYS
  REAL(prcn),ALLOCATABLE, DIMENSION(:,:) :: XC_PHYS
  REAL(prcn) :: numdens

  PUBLIC ::  post_mean_drag_poly, post_static_particles
  
  
CONTAINS

  SUBROUTINE post_static_particles(imis)
    IMPLICIT NONE
    INTEGER, Intent(in) :: imis
    INTEGER :: m, idim, imeasvol
    REAL(prcn) :: Lmeasmin, Lmeasmax, Lmeaslin, dmeasvol
    
    if(I_AM_NODE_ZERO) WRITE(*,'(A,I)') 'COMPUTING THE NEW TIME STEP FOR MIS',imis
    CALL compute_new_timestep(itrmax)
    CALL nonlinear(itrmax)
    
    CALL post_mean_drag(imis)
    
    numdens = REAL(nbody,prcn)/(DOML(1)*DOML(2)*DOML(3))

    if(ALLOCATED(XC_PHYS))then
       DEALLOCATE(XC_PHYS, RAD_PHYS)
    end if
       ALLOCATE(XC_PHYS(nbody,ndim), RAD_PHYS(nbody))
    !Convert Positions and radii to physical units
       Write(*,*)'nbody = ', nbody
    do m = 1, nbody
       do idim = 1, ndim
          XC_PHYS(m,idim) = (XC(m,idim)-one)*dx
       end do
       RAD_PHYS(m) = RADBDY(m)*dx
    end do
    Lmeasmin = 5.d0*dx
    Lmeasmax = DOML(2)

    dmeasvol = (Lmeasmax - Lmeasmin)/real(nmeasvols-1, prcn)

    do imeasvol = 1, nmeasvols
       Lmeaslin = Lmeasmin + (imeasvol-1) * dmeasvol
       Write(*,*)'Input measurement volume =', Lmeaslin
       CALL compute_num_vol_in_measvol(Lmeaslin, imeasvol, imis)
       LMEAS_VOLS(imeasvol) = Lmeaslin
       
    end do

    
  END SUBROUTINE post_static_particles


  SUBROUTINE compute_num_vol_in_measvol(Lmeaslin, imeasvol, imis)
    USE constants  
    IMPLICIT NONE
    INTEGER, Intent(in) :: imeasvol, imis
    REAL(prcn), Intent(in) :: Lmeaslin
    REAL(prcn) :: XW, XE, YS, YN, ZB, ZT, temp_pos(ndim), volume
    INTEGER :: m, number, idim
    INTEGER :: ifinew,ifinee, jfines, jfinen, kfineb,&
         & kfinet, ii, jj, kk, PIJK(ndim)
    LOGICAL :: XINMEASVOL, YINMEASVOL, ZINMEASVOL
    REAL(prcn) :: num_analy, alpha2, ssfm, homog_num, ssfm_vol, homog_vol, vol_anal

    Num_analy = numdens * (DOML(1)*DOML(2)*DOML(3))
    homog_num = numdens * (Lmeaslin**3.d0)
    
    !homog_vol = homog_num*pi*dia_phys**3.d0/6.d0
    homog_vol = homog_num*pi*(dbydx**3.d0)/6.d0
    
    XW = DOML(1)/two - Lmeaslin/two
    XE = DOML(1)/two + Lmeaslin/two
    
    YS = DOML(2)/two - Lmeaslin/two
    YN = DOML(2)/two + Lmeaslin/two

    ZB = DOML(3)/two - Lmeaslin/two
    ZT = DOML(3)/two + Lmeaslin/two

    number = 0
    do M = 1, nbody
       do idim = 1, ndim
          temp_pos(idim) = XC_PHYS(m,idim)
          if(temp_pos(idim).lt.zero)temp_pos(idim) = temp_pos(idim)+DOML(idim)
          if(temp_pos(idim).gt.DOML(idim))temp_pos(idim) = temp_pos(idim)-DOML(idim)
       end do
       
       XINMEASVOL = (temp_pos(1).ge.XW).and.(temp_pos(1).le.XE)
       YINMEASVOL = (temp_pos(2).ge.YS).and.(temp_pos(2).le.YN)
       ZINMEASVOL = (temp_pos(3).ge.ZB).and.(temp_pos(3).le.ZT)
       if(XINMEASVOL.AND.YINMEASVOL.AND.ZINMEASVOL)then
          number = number + 1
       else if(imeasvol.eq.nmeasvols)then
          Write(*,*)'Particles greater than DOML =', XINMEASVOL, YINMEASVOL, ZINMEASVOL
          Write(*,*)'ZY =', XC(m,3), ZB, ZT
       end if
    end do
    
    
    KFINEB = INT(ZB/dz) + 1
    KFINET = INT(ZT/dz) + 1
    JFINES = INT(YS/dy) + 1
    JFINEN = INT(YN/dy) + 1
    IFINEW = INT(XW/dx) + 1
    IFINEE = INT(XE/dx) + 1

    
    volume = zero
#if 0    
    do KK = KFINEB, KFINET
       do JJ = JFINES, JFINEN
          do II = IFINEW, IFINEE
             if((KK.le.mz).and.(JJ.le.my).and.(II.le.mx1))then
                if(.not.fluid_atijk(II,JJ,KK))then
                   !volume = volume + dx**3.d0
                   volume = volume + one**3.d0
                end if
             end if
          end do
       end do
    end do
#endif    
    
    NUMBER_IN_MIS(imis,imeasvol) =  real(number,prcn)/num_analy
    
    ssfm = (number**2.d0 - number)/(homog_num)**2.d0
    
    ssfm_vol = (volume**2.d0 - volume)/(homog_vol)**2.d0

#if 1
    if(imeasvol.eq.nmeasvols) then
       Write(*,*) ' Number in meas vol = ', number, num_analy
       Write(*,*) ' SSFM = ', ssfm
       Write(*,*) ' Volume in meas vol = ', volume, homog_vol
       Write(*,*) ' SSFM = ', ssfm_vol

    end if
#endif

    
    NUMBERSQ_IN_MIS(imis,imeasvol) = ssfm
    
    VOLUME_IN_MIS(imis,imeasvol) = volume/homog_vol

    VOLUMESQ_IN_MIS(imis,imeasvol) = ssfm_vol

  END SUBROUTINE compute_num_vol_in_measvol

  SUBROUTINE post_mean_drag(imis)
    IMPLICIT NONE
    INTEGER, Intent(in) :: imis
    INTEGER :: iphs

    do iphs = 1, nphases
       mis_norm_drag_spec(imis,iphs) = norm_drag_spec(iphs)
       mis_norm_drag_chem_spec(imis,iphs) = norm_drag_chem_spec(iphs)
    end do

  END SUBROUTINE post_mean_drag

  SUBROUTINE post_mean_drag_poly
    IMPLICIT NONE

    INTEGER :: normdragunit, iphs, timestepcount
    CHARACTER(LEN=80) :: FILENAME
    REAL(prcn) :: normdrag, normdragspec(nphases), tbytconv, tbytvis, tbytdiff, tmp_ferror_array(1:nerr_steps), mean_normdrag(nphases)
    LOGICAL :: FIRST_LINE, flowconverged

#if 0
    normdragunit = getnewunit(minunitno,maxunitno)
    
    FIRST_LINE = .TRUE.
    
    mean_normdrag = zero
    timestepcount = 0
    flowconverged = .FALSE.

    if(ALLOCATED(ferror_array)) DEALLOCATE(ferror_array)
    ALLOCATE(ferror_array(nerr_steps))

    do iphs = 1, nphases
       if(ASSOCIATED(phase_array(iphs)%ferror_array)) DEALLOCATE(phase_array(iphs)%ferror_array)
       ALLOCATE(phase_array(iphs)%ferror_array(nerr_steps))
    end do
    ferror_array(1:nerr_steps) = one
    do iphs = 1, nphases 
       phase_array(iphs)%ferror_array(1:nerr_steps) = one
    end do
    
    FILENAME = TRIM(RUN_NAME)//'_norm_drag.dat'
    OPEN(normdragunit,FILE=TRIM(FILENAME),status='old')
125 continue
    READ(normdragunit,*,END=225) tbytconv, tbytvis, tbytdiff, normdrag, normdragspec(1:nphases)
    !, usmean_act(1:ndim), (ufmean(idim)/ucharmod,idim=1,ndim)
    if(FIRST_LINE) then
       do iphs = 1, nphases
          phase_array(iphs)%ferror = ONE
       end do
       ferror = ONE
       FIRST_LINE = .FALSE.
    else
       do iphs = 1, nphases
          IF(normdragspec(iphs).gt.ZERO)then
             phase_array(iphs)%ferror = ABS(normdragspec(iphs) - phase_array(iphs)%fold)/normdragspec(iphs)
          ELSE
             phase_array(iphs)%ferror = ONE
          END IF
       end do
       if(normdrag.gt.zero)then
          ferror = ABS(normdrag - fold)/normdrag
       else
          ferror = ONE
       end if
    end if

    do iphs = 1, nphases
       phase_array(iphs)%fold = normdragspec(iphs)
    end do
    fold = normdrag

    tmp_ferror_array(1:nerr_steps) = ferror_array(1:nerr_steps)
    ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
    ferror_array(1) = ferror
    ferror_hist = SUM(ferror_array(1:nerr_steps))/real(nerr_steps,prcn)

    do iphs = 1, nphases
       tmp_ferror_array(1:nerr_steps) = phase_array(iphs)%ferror_array(1:nerr_steps)
       phase_array(iphs)%ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
       phase_array(iphs)%ferror_array(1) = phase_array(iphs)%ferror
!!$       PRINT*,'FERROR_A =', FERROR_ARRAY
       phase_array(iphs)%ferror_hist = SUM(phase_array(iphs)%ferror_array(1:nerr_steps))/real(nerr_steps,prcn)
    end do
    
    if(.not.flowconverged)then
       flowconverged = ferror_hist.lt.tol_ferror
#if 0       
       flowconverged = .TRUE.
       do iphs = 1, nphases
          flowconverged = flowconverged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
       end do
#endif
       if(flowconverged)then
          wRITE(*,*)'Ding Dong = ', timestepcount, ferror_hist, tol_ferror
          OPEN(unit = 1000, file=TRIM(RUN_NAME)//'_CONVERGED_NOW', status="replace")
          close(1000, status="keep")
       end if
    end if
    
    if(flowconverged) then
       do iphs = 1, nphases
          mean_normdrag(iphs) = mean_normdrag(iphs) + normdragspec(iphs)
       end do
       timestepcount = timestepcount + 1
       wRITE(*,*)' Normdrag = ', timestepcount, normdrag!phase_array(1:nphases)%ferror_hist, tol_ferror
    end if
    
    goto 125
    
225 continue
    Write(*,*)'FERROR HIST = ', ferror_hist, phase_array(1:nphases)%ferror_hist
    if(flowconverged)then
       wRITE(*,*)' Timestepcount = ', timestepcount, nerr_steps
       mean_normdrag(1:nphases) = mean_normdrag(1:nphases)/real(timestepcount,prcn)
       convmiscount = convmiscount + 1
       wRITE(*,*)' Normdrag = ', mean_normdrag(1:nphases)
       do iphs = 1, nphases
          mis_avg_normdrag(iphs) = mis_avg_normdrag(iphs) + mean_normdrag(iphs)
          stddev_normdrag(iphs) = stddev_normdrag(iphs) + mean_normdrag(iphs)**2.d0
       end do
    end if
#endif    
  END SUBROUTINE post_mean_drag_poly

END MODULE post_static
