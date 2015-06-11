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

MODULE lang_source_diss

  USE precision 
  USE constants 
  USE post_global_data
  USE general_funcs
  USE postproc_funcs
  USE randomno
  USE dependent_functions
  IMPLICIT NONE
#if 0

CONTAINS


  SUBROUTINE calc_source_diss
    IMPLICIT NONE  
    INTEGER :: funit,m, idim, unit1,ibin, jdim,nbody2, unit2, unit3,sunit, mark,i, iphs

    INTEGER, PARAMETER :: nhbins=30

    REAL(prcn) :: wt(nbody), ftemp(nbody), hist(nhbins), velmfcs(nbody,ndim), fluctv(nbody,ndim), forcemfcs(nbody,ndim),&
         & mean_vel(nphases,ndim), sdev_fluctv(nphases,ndim),dummy,&
         & accltn(nbody,ndim),flucta(nbody,ndim), mean_accltn(nphases,ndim),sdev_flucta(nphases,ndim) , mean_force(nphases,ndim), normf(nphases), mixmean

    REAL(prcn) :: fmin, fmax, wt1(nspec1), wt2(nspec2), u(nbody),ave,adev,sdev,var,skew,curt
    REAL(prcn) ::  gamma(nphases,ndim),source(nphases,ndim),diss(nphases,ndim), gamma_ls(nphases,ndim2), A(ndim,ndim),beta(nphases,ndim)
    Real(prcn) :: aivj(nphases,ndim,ndim), betgamchar(nphases),lchar(nphases),forchar(nphases),vfrac(nphases)
    REAL(prcn) :: beta_par(nphases),beta_per(nphases), gamma_par(nphases), gamma_per(nphases), source_par(nphases), source_per(nphases), diss_par(nphases), diss_per(nphases), gran_energy_out, meanacc(nphases)
    REAL(prcn) :: tkef,sdevu

!!$      DOUBLE PRECISION :: force_koch(nbody,3), accltn_koch(nbody,3), Re_p(3),flucta_koch(nbody,3), slope_koch, numer_koch,sdev_flucta_koch(3)

    INTEGER :: count_tmp
    CHARACTER*80 :: FILENAME
    CHARACTER*80 :: FILENAME1


    CALL screen_separator(80,'*')
    WRITE(*,*) 'READING THE FORMATTED FORCE DATA'

    funit = getnewunit(minunitno,maxunitno)
    !FILENAME = TRIM(RUN_NAME)//'_RESTART'
    !open(unit=funit,file=FILENAME,form="formatted"      &
    !      ,status="unknown")    
    !read(funit, *) count_tmp
    !close(funit, status='keep')
!!$    count_tmp = 1
!!$    WRITE(FILENAME1,'(I1)')count_tmp
!!$
!!$    FILENAME = TRIM(RUN_NAME)//'_FORCE_'//TRIM(FILENAME1)//'.dat'
!!$
!!$    open(unit=funit,file=FILENAME,form="formatted")
!!$    WRITE(*,*) 'NUMBER OF BODIES = ',nspec1, nspec2, nbody      
!!$    READ(funit,*)
!!$    READ(funit,*)
!!$
!!$    do m = 1, nbody
!!$       read(funit,*) force(m,1), force(m,2), force(m,3),force_insolid(m,1), force_insolid(m,2), force_insolid(m,3)
!!$       !read(funit,*) xc(m,1),xc(m,2), xc(m,3),force(m,1), force(m,2), force(m,3)
!!$    end do
!!$    close(funit, status='keep')
    !CALL compute_orthogonal_transform(A(1:ndim,1:ndim),uchar(1:ndim))
    
    !CALL transform_vector(A(1:ndim,1:ndim),force(1:nbody,1:ndim),forcemfcs(1:nbody,1:ndim),nbody)
    
    !CALL transform_vector(A(1:ndim,1:ndim),velbdy(1:nbody,1:ndim),velmfcs(1:nbody,1:ndim),nbody)
    
    ! Force computed from Koch's drag law
!!$    do m=1, nbody
!!$       Re_p(:) = (uchar(:)-velbdy(m,:))*dia_phys/(2.d0*vis)
!!$       do idim=1, 3
!!$       force_koch(m,idim) = (0.11519+0.180501*abs(Re_p(idim)))*(1-vol_frac1)*3*pi*vis*dia_phys*(uchar(idim)-velbdy(m,idim))! for vol frac = 0.2
!!$       enddo
!!$
!!$    end do

    !PRINT*,'FORCE =', sum(force(1:nspec1,1))/nspec1, sum(force(1:nspec1,2))/nspec1, sum(force(1:nspec1,3))/nspec1 
    beta = zero
    mean_force = zero
    normf = zero
    mixmean = zero
    vfrac = zero
    if(ibidisperse)then
       PRINT*,'DIAMETERS dia1 and dia2 = ', dia1, dia2
       PRINT*,'DIAMETER RATIO = ', dia2/DIA1
       lchar(1) = dia1
       lchar(2) = dia2
       vfrac(1) = vol_frac1
       vfrac(2) = vol_frac2
       do m=1,nspec1
          wt1(m) = 1.d0/real(nspec1, prcn)
          mean_force(1,:) = mean_force(1,:) + forcemfcs(m,:)*wt1(m)
          accltn(m,:) = forcemfcs(m,:)
!!$          accltn(m,:) = 6*forcemfcs(m,:)/(pi*dia1**3.d0)
          beta(1,:) = beta(1,:)+accltn(m,:)*wt1(m)
       enddo
       do m=1,nspec2
          wt2(m) = 1.d0/real(nspec2, prcn)
          mean_force(2,:) = mean_force(2,:) + forcemfcs(m,:)*wt2(m)
          accltn(nspec1+m,:) = forcemfcs(nspec1+m,:)
!!$          accltn(nspec1+m,:) = 6*forcemfcs(nspec1+m,:)/(pi*dia2**3.d0)
          beta(2,:) = beta(2,:)+accltn(m,:)*wt2(m)
       enddo
    else
       lchar(1) = dia_phys
       lchar(2) = one
       vfrac(1) = vol_frac1
       vfrac(2) = zero
       
       do m=1,nbody
          wt1(m) = 1.d0/real(nbody, prcn)
          mean_force(1,:) = mean_force(1,:) + force(m,:)*wt1(m)
          !accltn(m,:) = forcemfcs(m,:)
          accltn(m,:) = force(m,:)
!!$          accltn(m,:) = 6*forcemfcs(m,:)/(pi*dia_phys**3.d0)
!!$          accltn_koch(m,:) = 6*force_koch(m,:)/(pi*dia_phys**3.d0)
          !beta(1,:) = beta(1,:)+ accltn(m,:)*wt1(m)
       enddo
    endif
    
    do iphs = 1, nphases
       !beta_par(iphs) = beta(iphs,1)
       !beta_per(iphs) = half*(beta(iphs,2)+beta(iphs,3))
       forchar(iphs) = 3.d0*pi*vis*lchar(iphs)*ucharmod
       !forchar(iphs) = one
       normf(iphs) = DOT_PRODUCT(mean_force(iphs,1:ndim),mean_force(iphs,1:ndim))!DOT_PRODUCT(mean_force(iphs,1:ndim),uchar(1:ndim))
       normf(iphs) = DSQRT(normf(iphs))
       mixmean = mixmean + vfrac(iphs)*normf(iphs)
       
       ens_mean_force(iphs) = ens_mean_force(iphs) + normf(iphs)/forchar(iphs)
       ens_force_var(iphs) = ens_force_var(iphs) + (normf(iphs)/forchar(iphs))**2.d0
       
       !betgamchar(iphs) = 18.d0*vis/(lchar(iphs)**2.d0)
       !betgamchar(iphs) = one
       !beta_avg(iphs,1) = beta_avg(iphs,1) + beta_par(iphs)/betgamchar(iphs)
       !beta_avg(iphs,2) = beta_avg(iphs,2) + beta_per(iphs)/betgamchar(iphs)
       !beta_var(iphs,1) = beta_var(iphs,1) + (beta_par(iphs)/betgamchar(iphs))**2.d0
       !beta_var(iphs,2) = beta_var(iphs,2) + (beta_per(iphs)/betgamchar(iphs))**2.d0
    end do

    mixmean = mixmean/SUM(vfrac(1:nphases))
    mixmean = mixmean/(3.d0*pi*vis*ucharmod*dia_phys)
    mixmean_for =  mixmean_for + mixmean
    mixvar_for = mixvar_for + mixmean**2.d0
!!$    do idim=1,ndim
!!$       gamma(:,idim) = zero
!!$       source(:,idim) = zero
!!$       diss(:,idim) = zero
!!$    end do
    !gran_energy_out=zero

    if(.not.ibidisperse)then
       do idim=1, 3
          u(1:nbody) = accltn(1:nbody,idim)
          CALL moment1(4, nbody, nbody, wt1, u, mean_accltn(1,idim),adev,sdev_flucta(1,idim),var,skew,curt)
          flucta(1:nbody,idim) = (accltn(1:nbody,idim)-mean_accltn(1,idim))

          u(1:nbody) = velmfcs(1:nbody,idim)
          u(1:nbody) = velbdy(1:nbody,idim)
          CALL moment1(4, nbody, nbody, wt1, u, mean_vel(1,idim),adev,sdev_fluctv(1,idim),var,skew,curt)
          fluctv(1:nbody,idim) = (velmfcs(1:nbody,idim)-mean_vel(1,idim))
          fluctv(1:nbody,idim) = (velbdy(1:nbody,idim)-mean_vel(1,idim))
          !gran_energy_out= gran_energy_out + sdev_fluctv(1,idim)**2.d0
       end do

       PRINT*,'MEAN_ACCLTN = ', mean_accltn(1,1:ndim)   
       !       do i=1,nbody
       !         gran_energy_out = gran_energy_out + (velmfcs(i,1)&
       !               &-mean_vel(1,1))**2.d0 + (velmfcs(i,2)-mean_vel(1,2))&
       !             &**2.d0+ (velmfcs(i,3)-mean_vel(1,3))**2.d0
       !     end do
       ! gran_energy_out = (one/three)*gran_energy_out!/nbody
       ! PRINT*,'ReT out = ', sqrt(gran_energy_out)*lchar(1)/vis
    else

       do idim=1, 3
          u(1:nspec1) = accltn(1:nspec1,idim)
          CALL moment1(4, nspec1, nspec1, wt1, u, mean_accltn(1,idim),adev,sdev_flucta(1,idim),var,skew,curt)
          flucta(1:nspec1,idim) = (accltn(1:nspec1,idim)-mean_accltn(1,idim))

          u(1:nspec2) = accltn(nspec1+1:nspec1+nspec2,idim)
          CALL moment1(4, nspec2, nspec2, wt2, u, mean_accltn(2,idim),adev,sdev_flucta(2,idim),var,skew,curt)
          flucta(1:nspec1+nspec2, idim) = (accltn(1:nspec1+nspec2,idim)-mean_accltn(2,idim))
          
          u(1:nspec1) = velmfcs(1:nspec1,idim)
          CALL moment1(4, nspec1, nspec1, wt1, u, mean_vel(1,idim),adev,sdev_fluctv(1,idim),var,skew,curt)
          fluctv(1:nspec1,idim) = (velmfcs(1:nspec1,idim)-mean_vel(1,idim))
          
          u(1:nspec2) = velmfcs(nspec1+1:nspec1+nspec2,idim)
          CALL moment1(4, nspec2, nspec2, wt2, u, mean_vel(2,idim),adev,sdev_fluctv(2,idim),var,skew,curt)
          fluctv(nspec1+1:nspec1+nspec2,idim) = (velmfcs(nspec1+1:nspec1+nspec2,idim)- mean_vel(2,idim))
       end do
    end if
!!$    vprime_char = zero
    
!!$    do iphs = 1,nphases
!!$       meanacc(iphs) = DOT_PRODUCT(mean_accltn(iphs,1:ndim),mean_accltn(iphs,1:ndim))
!!$       meanacc(iphs)= sqrt(meanacc(iphs))
!!$       if(meanacc(iphs).gt.zero)then
!!$          do idim=1,ndim
!!$             ens_sdev_acc(iphs,idim) = ens_sdev_acc(iphs,idim) + sdev_flucta(iphs,idim)/meanacc(iphs)
!!$             sdev_acc_var(iphs,idim) = sdev_acc_var(iphs,idim) + (sdev_flucta(iphs,idim)/meanacc(iphs))**2.d0
!!$          end do
!!$       end if
!!$       vprime_char(iphs) = DOT_PRODUCT(sdev_fluctv(iphs,1:ndim),sdev_fluctv(iphs,1:ndim))
!!$       vprime_char(iphs) = sqrt(vprime_char(iphs)/three)
!!$       vprime_char(iphs) = one
!!$       !aprime_char(iphs) = betgamchar(iphs)*vprime_char(iphs)
!!$       
!!$       aprime_char(iphs) = DOT_PRODUCT(sdev_flucta(iphs,1:ndim),sdev_flucta(iphs,1:ndim))
!!$       aprime_char(iphs) = DSQRT(aprime_char(iphs))
!!$       aprime_char(iphs) = one
!!$    end do
       
!!$       OPEN(funit,FILE=TRIM(RUN_NAME)//"_Aivj_scatter.dat", status='replace')
!!$    do idim=1,ndim
!!$       write(funit,*)'ZONE'
!!$       do m=1,nspec1
!!$          write(funit,'(2(2x,f12.8),2x,(I5))') fluctv(m,idim)/sdev_fluctv(1,idim),flucta(m,idim)/sdev_flucta(1,idim),1
!!$       enddo
!!$       do m=nspec1+1,nbody
!!$          write(funit,'(2(2x,f12.8),2x,(I5))') fluctv(m,idim)/sdev_fluctv(2,idim),flucta(m,idim)/sdev_flucta(2,idim),-1
!!$       enddo
!!$    end do
!!$    close(funit,status='keep')
    ! Compute the Tensor <Aivj>
       aivj = zero
       CALL compute_aivj(aivj(1,1:ndim,1:ndim),force(1:nspec1,1:ndim),force(1:nspec1,1:ndim),nspec1)

!!$    if(nspec2.gt.0) then 
!!$       CALL compute_aivj(aivj(2,1:ndim,1:ndim),fluctv(nspec1+1:nbody,1:ndim),flucta(nspec1+1:nbody,1:ndim),nspec2)
!!$    end if
  
    do iphs = 1, 1
!       if(vprime_char(iphs).gt.zero)then
       do idim = 1, ndim
          do jdim = 1, ndim
             ens_aivj(iphs,idim,jdim) = ens_aivj(iphs,idim,jdim) + aivj(iphs,idim,jdim)/(sdev_flucta(iphs,idim)*sdev_flucta(iphs,jdim))   
             ens_aivj_var(iphs,idim,jdim) = ens_aivj_var(iphs,idim,jdim) + (aivj(iphs,idim,jdim)/(sdev_flucta(iphs,idim)*sdev_flucta(iphs,jdim)))**2.d0   
             
          end do
       end do
    end do
    
    

    ! Calculate the anisotropy in the tensor
    !CALL calc_anisotropy(aivj(1,1:ndim,1:ndim))
    !if(nspec2.gt.0)CALL calc_anisotropy(aivj(2,1:ndim,1:ndim))

    !CALL estimate_gamma_cartesian(fluctv(1:nspec1,1:ndim),flucta(1:nspec1,1:ndim),nspec1,gamma_ls(1,1:ndim2),ndim2)
    !if(nspec2.gt.0) CALL estimate_gamma_cartesian(fluctv(nspec1+1:nbody,1:ndim),flucta(nspec1+1:nbody,1:ndim),nspec2,gamma_ls(2,1:ndim2),ndim2)


    !CALL estimate_gamma_mfcs(fluctv(1:nspec1,1:ndim),flucta(1:nspec1,1:ndim),gamma(1,1:ndim),nspec1)

    !if(nspec2.gt.0) then 
    !   CALL estimate_gamma_mfcs(fluctv(nspec1+1:nbody,1:ndim),flucta(nspec1+1:nbody,1:ndim),gamma(2,1:ndim),nspec2)
    !end if
   
!!$    do iphs = 1, nphases
!!$       gamma_par(iphs) = gamma(iphs,1)
!!$       gamma_per(iphs) = half*(gamma(iphs,2)+gamma(iphs,3))
!!$       gamma_avg(iphs,1) = gamma_avg(iphs,1)+gamma_par(iphs)/betgamchar(iphs)
!!$       gamma_avg(iphs,2) = gamma_avg(iphs,2)+gamma_per(iphs)/betgamchar(iphs)
!!$
!!$       gamma_var(iphs,1) = gamma_var(iphs,1)+(gamma_par(iphs)/betgamchar(iphs))**2.d0
!!$       gamma_var(iphs,2) = gamma_var(iphs,2)+(gamma_per(iphs)/betgamchar(iphs))**2.d0
!!$    end do
!!$    PRINT*,'gamma_par = ', gamma_par
!!$    OPEN(funit, FILE=TRIM(RUN_NAME)//'_white-noise-pdf1.dat', form='formatted', status='replace')    
!!$
!!$    CALL compute_source_diss(fluctv(1:nspec1,1:ndim),flucta(1:nspec1,1:ndim),gamma(1,1:ndim),source(1,1:ndim),diss(1,1:ndim),nspec1,funit)
!!$    close(funit,status='keep')
!!$
!!$    if(nspec2.gt.0)then 
!!$       OPEN(funit, FILE=TRIM(RUN_NAME)//'_white-noise-pdf2.dat', form='formatted', status='replace')    
!!$       CALL compute_source_diss(fluctv(nspec1+1:nbody,1:ndim),flucta(nspec1+1:nbody,1:ndim),gamma(2,1:ndim),source(2,1:ndim),diss(2,1:ndim),nspec2,funit)
!!$       close(funit,status='keep')
!!$    end if
!!$
!!$    do iphs = 1, nphases
!!$       source_par(iphs) = source(iphs,1)
!!$       source_per(iphs) = half*(source(iphs,2)+source(iphs,3))
!!$       if(vprime_char(iphs).gt.zero)then
!!$          source_avg(iphs,1) = source_avg(iphs,1)+ source_par(iphs)/(aprime_char(iphs)*vprime_char(iphs))
!!$          source_avg(iphs,2) = source_avg(iphs,2)+ source_per(iphs)/(aprime_char(iphs)*vprime_char(iphs))
!!$
!!$          source_var(iphs,1) = source_var(iphs,1)+ (source_par(iphs)/(aprime_char(iphs)*vprime_char(iphs)))**2.d0
!!$          source_var(iphs,2) = source_var(iphs,2)+ (source_per(iphs)/(aprime_char(iphs)*vprime_char(iphs)))**2.d0
!!$       end if
!!$
!!$       diss_par(iphs) = diss(iphs,1)
!!$       diss_per(iphs) = half*(diss(iphs,2)+diss(iphs,3))
!!$       if(vprime_char(iphs).gt.zero)then
!!$          diss_avg(iphs,1) = diss_avg(iphs,1)+ diss_par(iphs)/(betgamchar(iphs)*vprime_char(iphs)**2.d0)
!!$          diss_avg(iphs,2) = diss_avg(iphs,2)+ diss_per(iphs)/(betgamchar(iphs)*vprime_char(iphs)**2.d0)
!!$
!!$          diss_var(iphs,1) = diss_var(iphs,1)+ (diss_par(iphs)/(betgamchar(iphs)*vprime_char(iphs)**2.d0))**2.d0
!!$          diss_var(iphs,2) = diss_var(iphs,2)+ (diss_per(iphs)/(betgamchar(iphs)*vprime_char(iphs)**2.d0))**2.d0
!!$       end if
!!$    end do


!!$    OPEN(funit,FILE=TRIM(RUN_NAME)//"_lang_model_params.dat", status='replace')
!!$    write(funit,*)'ZONE t=','"GAMMA"'
!!$    do idim=1,ndim
!!$       write(funit,'(2(2x,f12.8))')(gamma(i,idim),i=1,2)
!!$    end do
!!$
!!$    write(funit,*)'ZONE t=','"SOURCE"'
!!$    do idim=1,ndim
!!$       write(funit,'(2(2x,f12.8))')(source(i,idim),i=1,2)
!!$    end do
!!$
!!$    write(funit,*)'ZONE t=','"DISSIPATION"'
!!$    do idim=1,ndim
!!$       write(funit,'(2(2x,f12.8))')(diss(i,idim),i=1,2)
!!$    end do
!!$    close(funit,status='keep')
!!$    RETURN
!!$   CALL calc_fluid_tke(tkef,sdevu)
!!$   tke_fluid = tke_fluid + tkef!/(sdevu**two)
!!$   tke_fluid_var = tke_fluid_var + (tkef)!/(sdevu**two))**two
  end SUBROUTINE calc_source_diss

    


 SUBROUTINE calc_fluid_tke(tke,sdevu)
   USE dependent_functions
   USE nlmainarrays, velr=>ubcp
   IMPLICIT NONE
   REAL(prcn), Intent(out) :: tke,sdevu
   REAL(prcn) :: umean_int(ndim),ustdev(ndim)
   INTEGER :: i, j, k, idim
   CALL calc_velreal(velr)
   umean_int = zero
   ustdev = zero
   do i=1, mx1
       do k=1,mz
          do j=1,my
             if (fluid_atijk(i,j,k))then 
                do idim=1,ndim
                   umean_int(idim) = umean_int(idim) + velr(i,j,k,idim)
                   ustdev(idim) = ustdev(idim) + velr(i,j,k,idim)**two
                end do
             end if
             
          end do
       end do
    end do
    umean_int(:) = umean_int(:)/count_fluid
    ustdev(:) = DSQRT(ustdev(:)/count_fluid - umean_int(:)**two)
    !    sdevu = DOT_PRODUCT(ustdev(1:ndim),ustdev(1:ndim))
    sdevu = DOT_PRODUCT(umean_int(1:ndim),umean_int(1:ndim))
    sdevu = DSQRT(sdevu)
    tke = zero 
    do idim=1,ndim
      do k=1,mz
         do j=1,my
            do i=1,mx1
              if(fluid_atijk(i,j,k))then 
               tke = tke + (velr(i,j,k,idim)-umean_int(idim))**2.d0
             endif    
           end do
        end do
     end do
  end do
  tke = tke/(mx1*my*mz)
 END SUBROUTINE calc_fluid_tke 
 
 
 SUBROUTINE transform_vector(A,v,vhat,n)
   IMPLICIT NONE
   Integer, Intent(in) :: n
   REAL(prcn),Intent(in) :: v(n,ndim),A(ndim,ndim)
   REAL(prcn),Intent(inout) :: vhat(n,ndim)
   Integer :: m, idim
   REAL(prcn) :: vmfcs(ndim),vprime_mfcs(n,ndim),aprime_mfcs(n,ndim)
      
   do m = 1, n
      vmfcs(1:ndim) = v(m,1:ndim)
      vmfcs = matmul(A,vmfcs)
      vhat(m,1:ndim) = vmfcs(1:ndim)
   enddo
 END SUBROUTINE transform_vector


 SUBROUTINE estimate_gamma_mfcs(vprime_mfcs,aprime_mfcs,gamma,n)
   IMPLICIT NONE
   Integer, Intent(in) :: n
   REAL(prcn),Intent(inout) :: gamma(ndim)
   REAL(prcn),Intent(in) :: vprime_mfcs(n,ndim),aprime_mfcs(n,ndim)
   REAL(prcn) :: numer, denr,slope,adev,var,skew,curt
   Integer :: m, idim

   do idim=1,3
       numer=0
       denr=0
       do m=1,n
          numer=numer+ vprime_mfcs(m,idim)*aprime_mfcs(m,idim)
          denr= denr+ vprime_mfcs(m,idim)**2.d0
       end do

       slope = numer/denr
       gamma(idim) = abs(slope)             
       PRINT*,'gamma MFCS =',gamma(idim)
    end do

 END SUBROUTINE estimate_gamma_mfcs

 SUBROUTINE compute_source_diss(fluctv,flucta,gamma,source,diss,n,funit)
   IMPLICIT NONE
   Integer, Intent(in) :: n,funit
   REAL(prcn),Intent(in) :: fluctv(n,ndim),flucta(n,ndim),gamma(ndim)
   REAL(prcn),Intent(out) :: source(ndim),diss(ndim)
   Integer :: idim, m,nhbins,i
   parameter(nhbins=20)
   Real(prcn) :: ran_var(n), sigma_rv, mean_rv, stand_rv(n),u(n),wt(n), adev,var,skew, curt,ftemp(n),hist(nhbins), fmin, fmax

   diss = zero

   do idim = 1, ndim
      do m = 1, n
          diss(idim) = diss(idim) + fluctv(m,idim)**2.d0
       end do
       diss(idim) = two*gamma(idim)*diss(idim)/n
    end do

    do i = 1,n
       wt(i) = 1.d0/real(n,prcn)
    end do

    
    do idim = 1, ndim
       do m=1, n
          ran_var(m) = flucta(m,idim)+gamma(idim)*fluctv(m,idim)
       end do
       u(1:n) = ran_var(1:n)
       CALL moment1(4, n, n, wt, u, mean_rv,adev,sigma_rv,var,skew,curt)     
       source(idim) = sigma_rv**2.d0
       
       do m=1, n
          stand_rv(m) = (ran_var(m)-mean_rv)/sigma_rv
       end do
       ftemp(1:n) = stand_rv(1:n)
       CALL histogram(ftemp,wt,n,nhbins,fmin,fmax,hist)
       CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,funit,real(idim,prcn),1.d0)
       
    end do
    
    
  END SUBROUTINE compute_source_diss
 
 SUBROUTINE estimate_gamma_cartesian(fluctv,flucta,n,gamma_ls,ndim2)
   IMPLICIT NONE
   Integer, Intent(in) :: n,ndim2
   REAL(prcn),Intent(in) :: fluctv(n,ndim),flucta(n,ndim)
   REAL(prcn),Intent(inout) :: gamma_ls(ndim2)
   REAL(prcn), ALLOCATABLE, DIMENSION(:,:) :: vprime
   REAL(prcn), ALLOCATABLE, DIMENSION(:) :: aprime
   REAL(prcn)::  vprime_temp(ndim,ndim2)
   INTEGER :: idim,m,start, ending, row

   ALLOCATE(vprime(n*ndim, ndim2), aprime(n*ndim))
   gamma_ls(1:ndim2) = zero
   vprime = zero
   gamma_ls = zero
   aprime = zero
   vprime_temp = zero
   do  idim = 1, ndim
      start = (idim-1)*ndim + 1
      ending = start + ndim -1
      do m = 1, n
         row = (m-1)*ndim + idim
         vprime_temp(idim,start:ending) = fluctv(m,1:ndim)
         vprime(row,1:ndim2) = vprime_temp(idim,1:ndim2)
         aprime(row) = -flucta(m,idim)
      end do
   end do
   CALL estimate_gamma(vprime,aprime,gamma_ls,ndim*nbody, ndim2)
 END SUBROUTINE estimate_gamma_cartesian

 SUBROUTINE estimate_gamma(vprime,aprime,gamma,mv, nv)
   IMPLICIT NONE
!!$    This subroutine is used to estimate gamma by casting A' = -
   !!\gamma v' into a  linear least squares problem. The approximate
   !! set of equations to be solved are                             
   !!      V'\Gamma = -A' 
!!$  V' => mv X nv,  \Gamma => nv X 1,   A' => mv X 1
!!$ The equations are solved with the non negativity constraints in
   !! the last few rows of the \Gamma vector.

   INTEGER , INTENT(in) :: mv,nv
   INTEGER :: me,lw,liw,L,mdw,mode, row
   REAL(prcn), Intent(out) :: gamma(nv)
   REAL(prcn) :: prgopt(1),upval
   REAL(prcn),INTENT(inout) ::  vprime(mv,nv), aprime(mv)
   REAL(prcn)  :: rnorm
   REAL(prcn) :: vprime_temp(mv,nv), aprime_temp(mv),gamma_temp(nv)
   INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
   REAL(prcn), DIMENSION(:,:),ALLOCATABLE ::  E
   REAL(prcn), DIMENSION(:),ALLOCATABLE ::  F
   REAL(prcn), DIMENSION(:), ALLOCATABLE :: work
   REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: W


   me = 0 ! For the input to DWNNLS

!!$    Interchange rows and columns corresponding to \gamma_11, 
   !!\gamma_22 and \gamma_33 (since non negativity constraints are
   !! imposed on these values)
   vprime_temp = vprime
   do row = 1, mv
      vprime(row,1) = vprime_temp(row,nv-ndim+1) ! v'_i1 --> v'_i7
      vprime(row,nv-ndim+1) = vprime_temp(row,1) ! v'_i7 --> v'_i1
      vprime(row,1+ndim+1) = vprime_temp(row,nv-ndim+2) ! v'_i5 --> v'_i8
      vprime(row, nv-ndim+2) = vprime_temp(row,1+ndim+1) ! v'_i8 --> v'_i5
   end do
   aprime_temp = aprime
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
   W(me+1:me+mv,1:nv) = vprime(1:mv,1:nv)
   W(me+1:me+mv,nv+1) = aprime(1:mv)
   iwork(1) = lw
   iwork(2) = liw
   prgopt(1) = 1.
   !PRINT*,'MA = ', MV
   !PRINT*,'NA = ', NV
   !L = 3
   L=nv-ndim
   !PRINT*,'L = ', L
   CALL DWNNLS (W, MDW, ME, MV, NV, L, PRGOPT, gamma_temp, RNORM, MODE, &
        IWORK, WORK)
   gamma(1:nv) = gamma_temp(1:nv)
   !IF(mode.NE.0)   
   !PRINT*,'Mode=',mode

   gamma(1) = gamma_temp(nv-ndim+1)
   gamma(nv-ndim+1) = gamma_temp(1)
   gamma(1+ndim+1) = gamma_temp(nv-ndim+2)
   gamma(nv-ndim+2) = gamma_temp(1+ndim+1)
   !do row=1,nv
      !PRINT*,'gamma = ' ,gamma(row)
   !end do

 end SUBROUTINE estimate_gamma

 SUBROUTINE compute_orthogonal_transform(A,vel)
   IMPLICIT NONE

   Real(prcn),Intent(in) ::  vel(ndim)
   Real(prcn),Intent(inout) :: A(ndim,ndim)
   Real(prcn) ::  ct, st, velmag,theta,l1, l2, l3,rho
   Integer :: m
!!$
   velmag = vel(1)**2.d0 + vel(2)**2.d0 + vel(3)**2.d0
   velmag = sqrt(velmag)

   a(1,1) = vel(1)/velmag
   a(1,2) = vel(2)/velmag
   a(1,3) = vel(3)/velmag

   l1 = vel(1)/velmag
   l2 = vel(2)/velmag
   l3 = vel(3)/velmag

   PRINT*,'a(1,3)=', a(1,3)

   theta = pi/4.d0
   ct = cos(theta)
   PRINT*,'cos(theta) = ', ct
   st = sin(theta)
   
   rho = sqrt(l1**2.d0 + l2**2.d0)

   a(2,1) = (l1*l3*ct + l2*st)/rho
   a(2,2) = (l2*l3*ct - l1*st)/rho
   a(2,3) = -(rho**2.d0)*ct
   
   a(3,1) = l2*a(2,3)-l3*a(2,2)
   a(3,2) = l3*a(2,1)-l1*a(2,3)
   a(3,3) = l1*a(2,2)-l2*a(2,1)
 END SUBROUTINE compute_orthogonal_transform

 SUBROUTINE write_lang_param_output(countmis)
   IMPLICIT NONE
   INTEGER, Intent(in) :: countmis
   INTEGER :: idim,jdim,iphs, funit,j
   CHARACTER*80 :: FILENAME
   

   do idim=1,ndim
      do jdim = 1, ndim
         ens_aivj(1,idim,jdim) = ens_aivj(1,idim,jdim)/countmis
      end do
   enddo
   do idim=1,ndim
      do jdim = 1, ndim
         ens_aivj_var(1,idim,jdim) = ens_aivj_var(1,idim,jdim)/countmis - ens_aivj(1,idim,jdim)**2.d0
      end do
   enddo
!!$   ens_mean_force(:) = ens_mean_force(:)/countmis
!!$   ens_force_var(:) = ens_force_var(:)/countmis - (ens_mean_force(:))**2.d0
!!$   PRINT*,'en_force_var = ', ens_force_var
!!$   do idim=1,ndim1
!!$      beta_avg(:,idim) = beta_avg(:,idim)/countmis
!!$      beta_var(:,idim) = beta_var(:,idim)/countmis - beta_avg(:,idim)**2.d0
!!$      gamma_avg(:,idim) = gamma_avg(:,idim)/countmis
!!$      gamma_var(:,idim) = gamma_var(:,idim)/countmis - gamma_avg(:,idim)**2.d0
!!$      source_avg(:,idim) = source_avg(:,idim)/countmis
!!$      source_var(:,idim) = source_var(:,idim)/countmis - source_avg(:,idim)**2.d0
!!$      diss_avg(:,idim) = diss_avg(:,idim)/countmis
!!$      diss_var(:,idim) = diss_var(:,idim)/countmis - diss_avg(:,idim)**2.d0
!!$      ens_sdev_acc(:,idim) = ens_sdev_acc(:,idim)/countmis
!!$      sdev_acc_var(:,idim) = (sdev_acc_var(:,idim)/countmis) - ens_sdev_acc(:,idim)**2.d0
!!$   end do
   mixmean_for = mixmean_for/countmis
   mixvar_for = mixvar_for/countmis - mixmean_for**2.d0
   
!!$   tke_fluid = tke_fluid/countmis
!!$   tke_fluid_var = tke_fluid_var/countmis - tke_fluid**2.d0
   CALL screen_separator(80,'+')

!!$   FILENAME = TRIM(POST_RUNNAME)//'_beta_par.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(beta_avg(iphs,1),iphs=1,nphases),(confin*sqrt(beta_var(iphs,1)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_beta_per.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(beta_avg(iphs,2),iphs=1,nphases),(confin*sqrt(beta_var(iphs,2)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_gamma_par.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(gamma_avg(iphs,1),iphs=1,nphases),(confin*sqrt(gamma_var(iphs,1)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_gamma_per.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(gamma_avg(iphs,2),iphs=1,nphases),(confin*sqrt(gamma_var(iphs,2)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_source_par.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(source_avg(iphs,1),iphs=1,nphases),(confin*sqrt(source_var(iphs,1)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_source_per.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(source_avg(iphs,2),iphs=1,nphases),(confin*sqrt(source_var(iphs,2)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_diss_par.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(diss_avg(iphs,1),iphs=1,nphases),(confin*sqrt(diss_var(iphs,1)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_diss_per.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(diss_avg(iphs,2),iphs=1,nphases),(confin*sqrt(diss_var(iphs,2)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')

   FILENAME = TRIM(POST_RUNNAME)//'_mixmean_drag.dat'
   OPEN(funit,FILE=FILENAME,status='unknown')
   Write(funit,'(2(2x,g17.6))')mixmean_for, confin*sqrt(mixvar_for/countmis)
   Close(funit,status='keep')

!!$   FILENAME = TRIM(POST_RUNNAME)//'_mean_drag.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(4(2x,g17.6))')(ens_mean_force(iphs),iphs=1,nphases), (confin*sqrt(ens_force_var(iphs)/countmis),iphs=1,nphases)
!!$   Close(funit,status='keep')
!!$
!!$   FILENAME = TRIM(POST_RUNNAME)//'_acc_var.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   do idim=1,ndim
!!$      Write(funit,'(4(2x,g17.6))')(ens_sdev_acc(iphs,idim),iphs=1,nphases), (confin*sqrt(sdev_acc_var(iphs,idim)/countmis),iphs=1,nphases)
!!$   enddo
!!$   Close(funit,status='keep')

   FILENAME = TRIM(POST_RUNNAME)//'_fifj.dat'
   OPEN(funit,FILE=FILENAME,status='unknown')
   do idim=1,ndim
      Write(funit,'(6(2x,g17.6))')(ens_aivj(1,idim,j),j=1,ndim),(confin*ens_aivj_var(1,idim,j)/countmis, j=1,ndim)
   enddo
   Close(funit,status='keep')
     PRINT*,'FROM LANG_SOURCE_DISS.f90'
!!$    do idim = 1, ndim
!!$       PRINT*,(ENS_aivj(1,idim,j),j = 1, ndim)
!!$    end do
!!$   do iphs = 1,nphases
!!$      CALL calc_anisotropy(ens_aivj(iphs,1:ndim,1:ndim))
!!$   end do
!!$   FILENAME = TRIM(POST_RUNNAME)//'_fluid_tke.dat'
!!$   OPEN(funit,FILE=FILENAME,status='unknown')
!!$   Write(funit,'(2(2x,g17.6))')tke_fluid, confin*sqrt(tke_fluid_var/countmis)
!!$   Close(funit,status='keep')
   CALL screen_separator(80,'+')

 END SUBROUTINE write_lang_param_output
#endif   
end MODULE lang_source_diss

