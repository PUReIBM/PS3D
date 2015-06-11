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

MODULE langevin_model
#include "../FLO/ibm.h"
  USE precision 
  USE constants 
  USE post_global_data
  USE general_funcs
  USE postproc_funcs
  USE randomno
  IMPLICIT NONE

CONTAINS

  SUBROUTINE populate_current_realization(imis)
    IMPLICIT NONE
    INTEGER, Intent(in) :: imis
    INTEGER :: m, npart, idim, partstart, partend, iphs, jdim, funit
    REAL(prcn) :: mean_force(ndim), sigma_f(ndim), adev, var, skew,&
         & curt, wt1(nbody), u(nbody), ufvar(ndim), uiuj(ndim,ndim),&
         & mean_pres(ndim), mean_visc(ndim), sigma_pres(ndim),&
         & sigma_visc(ndim),force_koch(nbody,ndim),&
         & Re_p(ndim), Re_koch, avg_force_spec(nphases,ndim),&
         & avg_force(ndim), normdrag


    ! Populate current realization
    current_mis%mis = imis
    current_mis%nphases = nphases

    ALLOCATE(current_mis%phase_info(nphases))

    npart = 0

    do iphs = 1, nphases
       current_mis%phase_info(iphs)%dia = phase_array(iphs)%dia
       current_mis%phase_info(iphs)%npart = phase_array(iphs)%npart
       current_mis%phase_info(iphs)%volfrac = phase_array(iphs)%volfrac
       current_mis%phase_info(iphs)%volfracg = phase_array(iphs)%volfracg
       npart = npart + current_mis%phase_info(iphs)%npart
    end do

    current_mis%maxvolfrac = maxvolfrac
    ! ALLOCATE POINTER MEMORY
    
    IF (ASSOCIATED(current_mis%for)) THEN
       DEALLOCATE(current_mis%for)
       IF (npart.GT.0) then 
          ALLOCATE(current_mis%for(npart,ndim))
       end IF
    ELSE
       IF(npart.GT.0) ALLOCATE(current_mis%for(npart,ndim))
    ENDIF

    IF (ASSOCIATED(current_mis%vel)) THEN
       DEALLOCATE(current_mis%vel)
       IF (npart.GT.0) then 
          ALLOCATE(current_mis%vel(npart,ndim))
       end IF
    ELSE
       IF(npart.GT.0) ALLOCATE(current_mis%vel(npart,ndim))
    ENDIF

    
    do m = 1, npart
       current_mis%for(m,1:ndim) = force(m,1:ndim)
       current_mis%vel(m,1:ndim) = velbdy(m,1:ndim)

       do idim = 1, ndim
          if(ENSAVG%force_max(idim).lt.force(m,idim))then
             ENSAVG%force_max(idim) = force(m,idim)
          end if
          if(ENSAVG%force_min(idim).gt.force(m,idim))then
             ENSAVG%force_min(idim) = force(m,idim)
          end if

          if(ENSAVG%vel_max(idim).lt.velbdy(m,idim))then
             ENSAVG%vel_max(idim) = velbdy(m,idim)
          end if
          if(ENSAVG%vel_min(idim).gt.velbdy(m,idim))then
             ENSAVG%vel_min(idim) = velbdy(m,idim)
          end if
       end do
    end do
    
    do m = 1, nbody
       wt1(m) = 1/real(nbody,prcn)
    end do
    
    do idim=1, 3
       u(1:nbody) = force(1:nbody,idim)
       CALL moment1(4, nbody, nbody, wt1, u, mean_force(idim),adev,sigma_f(idim),var,skew,curt)

       u(1:nbody) = pres(1:nbody,idim)
       CALL moment1(4, nbody, nbody, wt1, u, mean_pres(idim),adev,sigma_pres(idim),var,skew,curt)
       u(1:nbody) = visc(1:nbody,idim)
       CALL moment1(4, nbody, nbody, wt1, u, mean_visc(idim),adev,sigma_visc(idim),var,skew,curt)
!!$       u(1:nbody) = velbdy(1:nbody,idim)
!!$       CALL moment1(4, nbody, nbody, wt1, u, mean_vel(1,idim),adev,sdev_v(1,idim),var,skew,curt)
    end do
    
    do idim = 1, ndim
       current_mis%usmean(idim) = usmean_des(idim)
       current_mis%ufmean(idim) = ufmean(idim)
       current_mis%mixmeanslip(idim) = (one-maxvolfrac)*(usmean(idim)-ufmean(idim))
       current_mis%mean_force(idim) = mean_force(idim)
       current_mis%mean_pres(idim) = mean_pres(idim)
       current_mis%mean_visc(idim) = mean_visc(idim)
       current_mis%sigma_f(idim) = sigma_f(idim)
       current_mis%sigma_pres(idim) = sigma_pres(idim)
       current_mis%sigma_visc(idim) = sigma_visc(idim)
    end do
    current_mis%gran_temp = gran_temp
    CALL calc_fluid_vel_statistics(ufvar(1:ndim), uiuj(1:ndim,1:ndim))

    do m=1, nbody
       Re_p(1:ndim) = (velbdy(m,1:ndim)-ufmean_des(1:ndim))*dia_phys&
            &/(2.d0*vis)
       Re_koch = DSQRT(DOT_PRODUCT(Re_p(1:ndim), Re_p(1:ndim)))
       do idim = 1, ndim
          force_koch(m,idim) = -(0.11519+0.180501*abs(Re_koch))*(1&
               &-maxvolfrac)*3*pi*vis*dia_phys*(velbdy(m,idim)-ufmean_des(idim))! for vol frac = 0.2
       enddo

    end do
    

    do idim = 1, ndim
       current_mis%ufvar(idim) = ufvar(idim)
       do jdim = 1, ndim
          current_mis%uiuj(idim,jdim) = uiuj(idim,jdim)
       end do
    end do
    !CALL calc_anisotropy(current_mis%uiuj(1:ndim,1:ndim), zi, eta)
    
    funit = getnewunit(minunitno,maxunitno)
    
    OPEN(unit=funit,FILE=TRIM(RUN_NAME)//'_fluid_anisotroy.dat', status&
         &='unknown')
    !write(funit,'(2(2x, g17.8))') zi, eta
    Close(funit, status='keep')
    
    OPEN(unit=funit,FILE=TRIM(RUN_NAME)//'_force_vel_scatter.dat', status&
         &='unknown')
    Do m = 1, nbody
       Write(funit,'(6(2x, g17.8))')(velbdy(m,idim), idim = 1, ndim),&
            & (force(m,idim), idim = 1, ndim)
    End Do
    Close(funit, status='keep')

    OPEN(unit=funit,FILE=TRIM(RUN_NAME)//'_force_koch_scatter.dat', status&
         &='unknown')
    
    Do m = 1, nbody
       Write(funit,'(6(2x, g17.8))')(velbdy(m,idim), idim = 1, ndim),&
            & (force_koch(m,idim), idim = 1, ndim)
    End Do

    Close(funit, status='keep')

    ALLOCATE(current_mis%next) ! Associated with unnamed storage
    current_mis => current_mis%next ! Finished working with current realization. Associate with the next realization
    
  END SUBROUTINE populate_current_realization
  
  SUBROUTINE compute_ensemble_avg_quantities
    IMPLICIT NONE
    INTEGER :: m, partstart, partend, npart(nphases),iphs&
         &,partcount(nphases),miscount,idim,jdim, funit, ibin
    REAL(prcn) :: fmean(ndim), fvar(ndim), pmean(ndim), pvar(ndim),&
         & grantemp, mixmeanslip(ndim), vismean(ndim), visvar(ndim)
    REAL(prcn) :: ufvar(ndim), sigma_ufvar(ndim)

    REAL(prcn) :: aivj(nphases,ndim,ndim),vivj(nphases,ndim,ndim),&
         & vfrac
    CHARACTER*80 :: FILENAME
    INTEGER, PARAMETER :: velbins = 100
    INTEGER :: count(velbins,ndim)
    REAL(prcn) :: vdiff(ndim), fhist(velbins, ndim),vhist(velbins,&
         & ndim), beta(velbins,ndim), fhistmag, vhistmag, normdrag, normdragvar, normtemp, scale
    
    fmean = zero
    fvar = zero
    
    pmean = zero
    pvar = zero
    
    vismean = zero
    visvar = zero
    
    ufvar = zero
    sigma_ufvar = zero
    
    mixmeanslip = zero
    grantemp = zero
    miscount = 0
    vfrac = zero
    
    current_mis => mis_data
    
    vdiff(1:ndim) = (ENSAVG%vel_max(1:ndim)-ENSAVG%vel_min(1:ndim))&
         &/real(velbins-1,prcn)
    do ibin = 1, velbins
       count(ibin,1:ndim) = 0
       fhist(ibin,1:ndim) = zero
       vhist(ibin,1:ndim) = zero
       beta(ibin,1:ndim) = zero
    end do
    normdrag = zero
    normdragvar = zero
    Do While(ASSOCIATED(current_mis%next))
       miscount = miscount + 1
       do idim = 1, ndim
          fmean(idim) = fmean(idim) + current_mis%mean_force(idim)
          fvar(idim) = fvar(idim) + current_mis%mean_force(idim)**2.d0
          
          pmean(idim) = pmean(idim) + current_mis%mean_pres(idim)
          pvar(idim) = pvar(idim) + current_mis%mean_pres(idim)**2.d0
          
          vismean(idim) = vismean(idim) + current_mis%mean_visc(idim)
          visvar(idim) = visvar(idim) + current_mis%mean_visc(idim)&
               &**2.d0
          
          mixmeanslip(idim) = mixmeanslip(idim) + current_mis&
               &%mixmeanslip(idim)
          grantemp = grantemp + current_mis%gran_temp
          
          ufvar(idim) = ufvar(idim) + current_mis%ufvar(idim)
          sigma_ufvar(idim) = sigma_ufvar(idim) + current_mis&
               &%ufvar(idim)**2.d0

             Write(*,*)'vdiff = ', idim, vdiff(idim)
          if(ABS(vdiff(idim)).gt.zero)then
            do m = 1, current_mis%phase_info(1)%npart
                ibin = (current_mis%vel(m,idim) - ENSAVG&
                       &%vel_min(idim))/vdiff(idim) + 1
                !Write(*,*)'ibin = ', ibin
                fhist(ibin,idim) = fhist(ibin,idim) + current_mis&
                     &%for(m,idim)
                vhist(ibin,idim) = vhist(ibin,idim) + current_mis&
                     &%vel(m,idim)
                count(ibin,idim) = count(ibin,idim) + 1
             end do
          end if
          normtemp = DSQRT(DOT_PRODUCT(current_mis%mean_force(1:ndim),current_mis%mean_force(1:ndim)))
          normdrag = normdrag + normtemp
          normdragvar = normdragvar + normtemp**two
       end do
       vfrac = vfrac + current_mis%maxvolfrac
       current_mis => current_mis%next
    end Do

    do ibin = 1, velbins
       do idim = 1, ndim
          if(count(ibin,idim).gt.0)then
             fhist(ibin,idim) = fhist(ibin,idim)/real(count(ibin&
                  &,idim),prcn)
             vhist(ibin,idim) = vhist(ibin,idim)/real(count(ibin&
                  &,idim),prcn)
             beta(ibin,idim) = -fhist(ibin,idim)/vhist(ibin,idim)
             
!!$          fhistmag = DOT_PRODUCT(fhist(ibin,1:ndim), vhist(ibin,&
!!$               & 1:ndim))
!!$          vhistmag = DOT_PRODUCT(vhist(ibin,1:ndim), vhist(ibin,&
!!$               & 1:ndim))
          end if
       end do
    end do

    funit = getnewunit(minunitno, maxunitno)    
    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_beta.dat', status =&
         & 'unknown')
    do ibin = 1, velbins
       Write(funit,'(6(2x,g17.8))')(vhist(ibin,idim), idim=1,ndim)&
            &,(beta(ibin,idim), idim=1,ndim)
    end do
    
    close(funit, status='keep')

    do idim = 1, ndim
       fmean(idim) = fmean(idim)/real(miscount,prcn)
       fvar(idim) = fvar(idim)/real(miscount,prcn) - fmean(idim)**2.d0
       
       pmean(idim) = pmean(idim)/real(miscount,prcn)
       pvar(idim) = pvar(idim)/real(miscount,prcn) - pmean(idim)**2.d0
       
       vismean(idim) = vismean(idim)/real(miscount,prcn)
       visvar(idim) = visvar(idim)/real(miscount,prcn) - vismean(idim)**2.d0
       
       ufvar(idim) = ufvar(idim)/real(miscount,prcn)
       sigma_ufvar(idim) = sigma_ufvar(idim)/real(miscount,prcn) - ufvar(idim)**2.d0
       
       mixmeanslip(idim) = mixmeanslip(idim)/real(miscount,prcn)
    end do
    normdrag = normdrag/real(miscount,prcn)
    normdragvar = normdragvar/real(miscount,prcn) - normdrag**2.d0
    
    vfrac = vfrac/real(miscount,prcn)
    grantemp = grantemp/real(miscount,prcn)    
    
    ENSAVG%nmis = miscount
    do idim = 1, ndim
       ENSAVG%mean_force(idim) = fmean(idim)
       ENSAVG%sigma_f(idim) = DSQRT(fvar(idim)) 
       
       ENSAVG%mean_pres(idim) = pmean(idim)
       ENSAVG%sigma_pres(idim) = DSQRT(pvar(idim)) 
       
       ENSAVG%mean_visc(idim) = vismean(idim)
       ENSAVG%sigma_visc(idim) = DSQRT(visvar(idim)) 
       
       ENSAVG%ufvar(idim) = ufvar(idim)
       ENSAVG%sigma_ufvar(idim) = DSQRT(sigma_ufvar(idim)) 
       
       ENSAVG%meanslip(idim) = mixmeanslip(idim)
    end do
    ENSAVG%grantemp = grantemp
    ENSAVG%volfrac = vfrac
    ENSAVG%meanslipmod = DSQRT(DOT_PRODUCT(mixmeanslip(1:ndim)&
         &,mixmeanslip(1:ndim)))
    ENSAVG%norm_drag = DSQRT(DOT_PRODUCT(fmean(1:ndim)&
         &,fmean(1:ndim)))
    
    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_normdragmod.dat',&
         & status = 'unknown')
    scale = 3.d0*pi*vis*ENSAVG%meanslipmod*dia_phys
    Write(funit,'(2(2x,g17.8))')(normdrag/scale), confin&
             &*normdragvar/(scale*real(ENSAVG&
            &%nmis,prcn))
    CLOSE(funit, status='keep')
  END SUBROUTINE compute_ensemble_avg_quantities

  SUBROUTINE Write_Ensemble_average_Data
    IMPLICIT NONE
    INTEGER :: funit, idim
    REAL(prcn) :: scale

    funit = getnewunit(minunitno, maxunitno)

    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_drag_unscaled.dat', status = 'unknown')
    Write(funit,'(18(2x,g17.8))')(ENSAVG%mean_force(idim), idim = 1,&
         & ndim),(ENSAVG%mean_pres(idim), idim = 1, ndim), (ENSAVG&
         &%mean_visc(idim), idim = 1, ndim), (confin*ENSAVG& 
         &%sigma_f(idim)/real(ENSAVG%nmis,prcn), idim = 1, ndim),&
         & (confin*ENSAVG& 
         &%sigma_pres(idim)/real(ENSAVG%nmis,prcn), idim = 1, ndim)&
         &,(confin*ENSAVG& 
         &%sigma_visc(idim)/real(ENSAVG%nmis,prcn), idim = 1, ndim)
    CLOSE(funit, status='keep')

    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_drag_scaled_slip.dat', status = 'unknown')
    
    scale = 3.d0*pi*vis*ENSAVG%meanslipmod*dia_phys
    Write(funit,'(18(2x,g17.8))')(ENSAVG%mean_force(idim)/scale, idim = 1,&
         & ndim),(ENSAVG%mean_pres(idim)/scale, idim = 1, ndim), (ENSAVG&
         &%mean_visc(idim)/scale, idim = 1, ndim), (confin*ENSAVG& 
         &%sigma_f(idim)/(scale*real(ENSAVG%nmis,prcn)), idim = 1, ndim),&
         & (confin*ENSAVG& 
         &%sigma_pres(idim)/(scale*real(ENSAVG%nmis,prcn)), idim = 1, ndim)&
         &,(confin*ENSAVG& 
         &%sigma_visc(idim)/(scale*real(ENSAVG%nmis,prcn)), idim = 1, ndim)
    
    CLOSE(funit, status='keep')
    
    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_drag_scaled_temp.dat', status = 'unknown')

    scale = 3.d0*pi*vis*DSQRT(ENSAVG%grantemp)*dia_phys
    Write(funit,'(18(2x,g17.8))')(ENSAVG%mean_force(idim)/scale, idim&
         & = 1,&
         & ndim),(ENSAVG%mean_pres(idim)/scale, idim = 1, ndim), (ENSAVG&
         &%mean_visc(idim)/scale, idim = 1, ndim), (confin*ENSAVG& 
         &%sigma_f(idim)/(scale*real(ENSAVG%nmis,prcn)), idim = 1, ndim),&
         & (confin*ENSAVG& 
         &%sigma_pres(idim)/(scale*real(ENSAVG%nmis,prcn)), idim = 1, ndim)&
         &,(confin*ENSAVG& 
         &%sigma_visc(idim)/(scale*real(ENSAVG%nmis,prcn)), idim = 1,&
         & ndim)
    
    CLOSE(funit, status='keep')

    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_ufvar.dat', status = 'unknown')
    
    scale = DSQRT(ENSAVG%grantemp)
    Write(funit,'(6(2x,g17.8))')(ENSAVG%ufvar(idim)/scale, idim&
         & = 1,ndim), (confin*ENSAVG%sigma_ufvar(idim)/(scale*real(ENSAVG&
         &%nmis,prcn)),idim = 1, ndim)
    
    CLOSE(funit, status='keep')

    OPEN(unit=funit, FILE = TRIM(POST_RUNNAME)//'_norm_drag.dat',&
         & status = 'unknown')
    scale = 3.d0*pi*vis*ENSAVG%meanslipmod*dia_phys
    Write(funit,'(6(2x,g17.8))')(ENSAVG%norm_drag/scale)!, (confin&
    !         &*ENSAVG%normdrag_var(idim)/(scale*real(ENSAVG&
    !        &%nmis,prcn)))
    CLOSE(funit, status='keep')

  END SUBROUTINE Write_Ensemble_average_Data

  SUBROUTINE calc_fluid_vel_statistics(ufvar,uiuj)
    USE dependent_functions
    USE nlmainarrays, velr=>ubcp

    IMPLICIT NONE
    REAL(prcn), Intent(out) :: ufvar(ndim), uiuj(ndim,ndim)
    REAL(prcn),DIMENSION(:),ALLOCATABLE :: fluid_vel, wt1
    Integer :: i,j,k,idim,flcount,funit, nhbins,m, jdim
    parameter(nhbins=100)
    REAL(prcn) :: fmin,fmax,hist(nhbins),usdev(ndim),adev,var,skew,curt
    
    
    ALLOCATE(wt1(count_fluid))
    
    do i = 1,count_fluid
       wt1(i) = 1.d0/real(count_fluid,prcn)
    end do
    
    funit = getnewunit(minunitno,maxunitno)

    OPEN(funit,FILE=TRIM(RUN_NAME)//'_fluid_vel_pdfs.dat', status='unknown')
    
    do idim = 1, ndim
       if(.not.ALLOCATED(fluid_vel))then
          ALLOCATE(fluid_vel(count_fluid))
       end if
       flcount = 1
       ufvar(idim) = zero
       do jdim = 1, ndim
          uiuj(idim,jdim) = zero
       end do
       
       do k=1,mz
          do j = 1, my
             do i = 1, mx1
                if(fluid_atijk(i,j,k))then
                   fluid_vel(flcount) = velr(i,j,k,idim) - ufmean(idim)
                   ufvar(idim) = ufvar(idim) + fluid_vel(flcount)**2.d0
                   do jdim = 1, ndim
                      uiuj(idim,jdim) = uiuj(idim,jdim) + (velr(i,j,k,idim)-ufmean(idim))*(velr(i,j,k,jdim)-ufmean(jdim))
                   end do
                   flcount = flcount + 1
                end if
             end do
          end do
       end do

       ufvar(idim) = ufvar(idim)/real(flcount,prcn)
       
       CALL histogram(fluid_vel,wt1,count_fluid,nhbins,fmin,fmax,hist)
       CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,funit,real(idim,prcn),1.d0)
    end do
    do idim = 1, ndim
       do jdim = 1, ndim
          uiuj(idim,jdim) = uiuj(idim,jdim)/real(count_fluid, prcn)
       end do
    end do

    DEALLOCATE(wt1)
          
  end SUBROUTINE calc_fluid_vel_statistics



#if 0
  SUBROUTINE estimate_model_coeffs_method1
    IMPLICIT NONE
    
    REAL(prcn) :: beta(nphases),fmean(nphases,ndim),mean_slip(nphases,ndim),slipmod(nphases),vprime_temp(ndim)
    INTEGER :: iphs,n,partstart,partend,row, start, ending,funit,nconst,jdim,idim,col,partcount,iiphs,dim1
    REAL(prcn), ALLOCATABLE, DIMENSION(:,:) :: vprime
    REAL(prcn), ALLOCATABLE, DIMENSION(:) :: aprime
    REAL(prcn) :: gamma_ls(ndim2),Rij(ndim,ndim),aivj(ndim,ndim)
    funit = getnewunit(minunitno,maxunitno)
    OPEN(unit=funit,FILE=TRIM(POST_RUNNAME)//'_LANG_COEFS_METHOD1.dat',form='formatted')
    
    
    beta = zero
    fmean(1:nphases,1:ndim) = ENSAVG%fmean(1:nphases,1:ndim)
    mean_slip(1:nphases,1:ndim) = ENSAVG%W(1:nphases,1:ndim)
    
    ! Estimate beta
    do iphs = 1,nphases
       beta(iphs) = DOT_PRODUCT(fmean(iphs,1:ndim),mean_slip(iphs,1:ndim))
       slipmod(iphs) = DOT_PRODUCT(mean_slip(iphs,1:ndim),mean_slip(iphs,1:ndim))
       if(slipmod(iphs).gt.zero)then
          beta(iphs) = beta(iphs)/slipmod(iphs)
       end if
    end do
    MODEL1%beta(1:nphases) = beta(1:nphases)
    
    do iphs = 1,nphases
       N = ENSAVG%nsamples(iphs)
       gamma_ls(1:ndim2) = zero
       do idim =1 ,ndim
          do jdim = 1, ndim
             Rij(idim,jdim) = ENSAVG%Rij(iphs,idim,jdim)
             aivj(idim,jdim) = ENSAVG%aivj(iphs,idim,jdim)
          end do
       end do

       if(N.gt.0)then
          ALLOCATE(vprime(n*ndim, ndim2), aprime(n*ndim))
          vprime = zero
          aprime = zero
          vprime_temp = zero
          row = 0
          
          current_mis => mis_data
          Do While(ASSOCIATED(current_mis%next))
             partstart = 1
             do iiphs = 1,iphs-1
                partstart = partstart + current_mis%npart(iiphs)
             end do
             partend = partstart+current_mis%npart(iphs)-1
             do partcount = partstart,partend
                vprime_temp(1:ndim) = current_mis%vel(partcount,1:ndim)
                do  idim = 1, ndim
                   row = row + 1
                   start = (idim-1)*ndim + 1
                   ending = start + ndim -1
                   vprime(row,start:ending) = vprime_temp(1:ndim)
                   aprime(row) = -current_mis%for(partcount,idim)
                end do
             end do
             current_mis => current_mis%next
          end Do
          nconst = 0
          CALL least_squares_fit(vprime,aprime,gamma_ls,N*ndim,ndim2,nconst)
          col = 0
          do idim = 1, ndim
             do jdim = 1, ndim
                col = col + 1
                MODEL1%gamma(iphs,idim,jdim) = gamma_ls(col)
             end do
          end do
          DEALLOCATE(vprime, aprime)
       end if
       do idim = 1,ndim
          do jdim = 1,ndim
             MODEL1%Bij(iphs,idim,jdim) = zero
          end do
       end do
       do idim = 1,ndim
          do jdim = idim,ndim
             do dim1 = 1,ndim
                MODEL1%Bij(iphs,idim,jdim) = MODEL1%Bij(iphs,idim,jdim) + MODEL1%gamma(iphs,idim,dim1)*Rij(jdim,dim1)
             end do
             do dim1 = 1,ndim
                MODEL1%Bij(iphs,idim,jdim) = MODEL1%Bij(iphs,idim,jdim) + MODEL1%gamma(iphs,jdim,dim1)*Rij(idim,dim1)
             end do
             MODEL1%Bij(iphs,idim,jdim) = MODEL1%Bij(iphs,idim,jdim) + aivj(idim,jdim)+aivj(jdim,idim)
             MODEL1%Bij(iphs,jdim,idim) = MODEL1%Bij(iphs,idim,jdim)
          end do
       end do
       
       Write(funit,*)'ZONE T = "SPECIES', iphs, ' "'
       do idim = 1, ndim
          Write(funit,'(6(2x,g17.6))')(MODEL1%gamma(iphs,idim,jdim),jdim=1,ndim),(MODEL1%Bij(iphs,idim,jdim),jdim=1,ndim)
       end do
    end do !NPHASES
  END SUBROUTINE estimate_model_coeffs_method1
  
  SUBROUTINE estimate_model_coeffs_method2
    IMPLICIT NONE
    REAL(prcn) :: beta(nphases),fmean(nphases,ndim),mean_slip(nphases,ndim),slipmod(nphases)
    REAL(prcn) :: Rij(nphases,ndim,ndim), aivj(nphases,ndim,ndim)
    REAL(prcn) :: COEF(6,15), sol(15), rhs(15)
    INTEGER :: iphs, idim, jdim, row, col, dim1, dim2, bijcol,nconst,funit
    
    funit = getnewunit(minunitno,maxunitno)
    OPEN(unit=funit,FILE=TRIM(POST_RUNNAME)//'_LANG_COEFS_METHOD2.dat',form='formatted')
    beta = zero
    fmean(1:nphases,1:ndim) = ENSAVG%fmean(1:nphases,1:ndim)
    mean_slip(1:nphases,1:ndim) = ENSAVG%W(1:nphases,1:ndim)
    
    ! Estimate beta
    do iphs = 1,nphases
       beta(iphs) = DOT_PRODUCT(fmean(iphs,1:ndim),mean_slip(iphs,1:ndim))
       slipmod(iphs) = DOT_PRODUCT(mean_slip(iphs,1:ndim),mean_slip(iphs,1:ndim))
       if(slipmod(iphs).gt.zero)then
          beta(iphs) = beta(iphs)/slipmod(iphs)
       end if
    end do
    MODEL2%beta(1:nphases) = beta(1:nphases)

    ! Estimate gamma_{ij} and B_{ij}
    do iphs = 1, ndim
       do idim =1 ,ndim
          do jdim = 1, ndim
             Rij(iphs,idim,jdim) = ENSAVG%Rij(iphs,idim,jdim)
             aivj(iphs,idim,jdim) = ENSAVG%aivj(iphs,idim,jdim)
          end do
       end do
    end do

    do iphs = 1, nphases
       ! Fill the first row of the coefficient matrix. 
       col = 0
       RHS = zero
       do idim = 1,ndim
          do jdim = 1,ndim
             col = col + 1 
             !COEF(1,col) = -two*Rij(iphs,idim,jdim)
             COEF(1,col) = two*Rij(iphs,idim,jdim)
          end do
          RHS(1) = RHS(1) + two * aivj(iphs,idim,idim)
       end do
       do idim = 1,ndim
          col = col + 1
          COEF(1,col) = zero
       end do
       do idim = 1,ndim
          col = col + 1
          COEF(1,col) = one
       end do

       ! Fill the subsequent rows of COEF.
       row = 0
       do idim = 1, ndim!ndim-1
          do jdim = idim, ndim
             row = row + 1
             col = 0
             COEF(row,1:15) = zero
             RHS(row) = zero
           !  PRINT*,'row = ', row
           !  if((idim.eq.jdim))then
           !     col = 0
           !     do dim1 = 1, ndim
           !        do dim2 = 1, ndim
           !           col = col + 1
           !           !COEF(row,col) = 2.d0/3.d0*Rij(iphs,dim1,dim2)
           !           COEF(row,col) = -2.d0/3.d0*Rij(iphs,dim1,dim2)
           !        enddo
           !     enddo
           !     do dim1 = 1,ndim
           !        col = (idim-1)*ndim + dim1
           !        !COEF(row,col) = COEF(row,col)-two*Rij(iphs,idim,dim1)
           !        COEF(row,col) = COEF(row,col)+two*Rij(iphs,idim,dim1)
           !        RHS(row) = RHS(row)-2.d0/3.d0*aivj(iphs,dim1,dim1)
           !     enddo
           !     col = ndim2 + ndim
           !     do dim1 = 1, ndim
           !        col = col + 1
           !        COEF(row,col) = -one/three
           !     enddo
           !     RHS(row) = RHS(row) +  two*aivj(iphs,idim,idim)
           !  else
                do dim1 = 1,ndim
                   col = (idim-1)*ndim + dim1
                   PRINT*,' col = ', col, jdim, dim1
                   !COEF(row,col) = -Rij(iphs,jdim,dim1)
                   COEF(row,col) = COEF(row,col)+Rij(iphs,jdim,dim1)
                end do
                do dim1 = 1,ndim
                   col = (jdim-1)*ndim + dim1
                   PRINT*,' col = ', col, idim, dim1
                   !COEF(row,col) = -Rij(iphs,idim,dim1)
                   COEF(row,col) = COEF(row,col)+Rij(iphs,idim,dim1)
                end do
                RHS(row) = aivj(iphs,idim,jdim)+aivj(iphs,jdim,idim)
             !end if
             bijcol = (idim-1)*ndim + jdim
             if(bijcol.eq.2.or.bijcol.eq.4)then
                bijcol = 1
             elseif(bijcol.eq.3.or.bijcol.eq.7)then
                bijcol = 2
             elseif(bijcol.eq.6.or.bijcol.eq.8)then
                bijcol = 3
             else
                bijcol = idim+ndim
             end if
             col = ndim*ndim + bijcol
             PRINT*,' col = ', col, bijcol
             COEF(row,col) = COEF(row,col) + one
          end do
       end do
       nconst = ndim
       CALL least_squares_fit(COEF,RHS,SOL,6,15,nconst)
       col = 0
       do idim = 1, ndim
          do jdim = 1, ndim
             col = col + 1
             MODEL2%gamma(iphs,idim,jdim) = SOL(col)
          end do
       end do
       do idim = 1,ndim
          do jdim = idim+1,ndim
             col = col+1
             MODEL2%Bij(iphs,idim,jdim) = SOL(col)
             MODEL2%Bij(iphs,jdim,idim) = MODEL2%Bij(iphs,idim,jdim)
          end do
       end do
       do idim =1 ,ndim
          do jdim = 1,ndim
             if(idim.eq.jdim)then
                col = col + 1
                MODEL2%Bij(iphs,idim,jdim) = SOL(col)
             end if
          end do
       end do
       do idim=10,15
          PRINT*,'SOL =', SOL(idim)
       end do
       Write(funit,*)'ZONE T = "SPECIES', iphs, ' "'
       do idim = 1, ndim
          Write(funit,'(6(2x,g17.6))')(MODEL2%gamma(iphs,idim,jdim),jdim=1,ndim),(MODEL2%Bij(iphs,idim,jdim),jdim=1,ndim)
       end do
    end do ! NPHASES
    CLOSE(funit,status='keep')
  END SUBROUTINE estimate_model_coeffs_method2

 SUBROUTINE least_squares_fit(COEF,RHS,SOL,nrow, ncol,nconst)
   IMPLICIT NONE
!!$    This subroutine is used to estimate gamma by casting A' = -
   !!\gamma v' into a  linear least squares problem. The approximate
   !! set of equations to be solved are                             
   !!      V'\Gamma = -A' 
!!$  COEF => nrow X ncol,  SOL => ncol X 1,   RHS => ncol X 1
!!$ The equations are solved with the non negativity constraints in
   !! the last few rows of the SOL vector. The input matrices should
!!$ be in such a way that the last rows should correspond to the
!!$   variables on which the non-negativity constraints are imposed

   INTEGER , INTENT(in) :: nrow,ncol,nconst
   INTEGER :: me,lw,liw,L,mdw,mode, row,mv,nv
   REAL(prcn), Intent(out) :: SOL(ncol)
   REAL(prcn) :: prgopt(1),upval
   REAL(prcn),INTENT(inout) ::  COEF(nrow,ncol), RHS(ncol)
   REAL(prcn)  :: rnorm, sol_temp(ncol)
   INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
   REAL(prcn), DIMENSION(:,:),ALLOCATABLE ::  E
   REAL(prcn), DIMENSION(:),ALLOCATABLE ::  F
   REAL(prcn), DIMENSION(:), ALLOCATABLE :: work
   REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: W
   
   
   
!!$    Interchange rows and columns corresponding to \gamma_11, 
   !!\gamma_22 and \gamma_33 (since non negativity constraints are
   !! imposed on these values)
!!$   vprime_temp = vprime
!!$   do row = 1, mv
!!$      vprime(row,1) = vprime_temp(row,nv-ndim+1) ! v'_i1 --> v'_i7
!!$      vprime(row,nv-ndim+1) = vprime_temp(row,1) ! v'_i7 --> v'_i1
!!$      vprime(row,1+ndim+1) = vprime_temp(row,nv-ndim+2) ! v'_i5 --> v'_i8
!!$      vprime(row, nv-ndim+2) = vprime_temp(row,1+ndim+1) ! v'_i8 --> v'_i5
!!$   end do
!!$   aprime_temp = aprime

   me = 0 ! For the input to DWNNLS
   mv = nrow
   nv = ncol
   mdw  = me+mv
   !k = max(ma+mg,nmat)
   lw = me+mv+5*nv
   liw = me+mv+nv
   IF(.NOT.ALLOCATED(iwork)) ALLOCATE(iwork(liw))
   IF(.NOT.ALLOCATED(work)) ALLOCATE(work(mv+5*nv))
   IF(.NOT.ALLOCATED(W)) ALLOCATE(W(mdw,nv+1))
   IF(me.GT.0) THEN 
      W(1:me,1:nv) = E(1:me,1:nv)
      W(1:me,nv+1) = F(1:me)
   ENDIF
   W(me+1:me+mv,1:nv) = COEF(1:mv,1:nv)
   W(me+1:me+mv,nv+1) = RHS(1:mv)
   iwork(1) = lw
   iwork(2) = liw
   prgopt(1) = 1.
   !PRINT*,'MA = ', MV
   !PRINT*,'NA = ', NV
   !L = 3
   L = ncol-nconst
!   L=nv-ndim
   !PRINT*,'L = ', L
   CALL DWNNLS (W, MDW, ME, MV, NV, L, PRGOPT, sol_temp, RNORM, MODE, &
        IWORK, WORK)
   SOL(1:nv) = sol_temp(1:nv)
   !IF(mode.NE.0)   
   !PRINT*,'Mode=',mode

!!$   gamma(1) = gamma_temp(nv-ndim+1)
!!$   gamma(nv-ndim+1) = gamma_temp(1)
!!$   gamma(1+ndim+1) = gamma_temp(nv-ndim+2)
!!$   gamma(nv-ndim+2) = gamma_temp(1+ndim+1)
   !do row=1,nv
      !PRINT*,'gamma = ' ,gamma(row)
   !end do
   
 end SUBROUTINE least_squares_fit

#endif
END MODULE langevin_model
