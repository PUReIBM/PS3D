!ROUTINE TO CALCULATE THE IMMERSED BOUNDARY FORCE FOR THE SCALAR FIELD
!Author: RAHUL GARG
MODULE bc_scalar
#include "../FLO/ibm.h"
  USE precision 
  USE constants 
  USE fftw_interface
  USE scalar_data

  Use nlmainarrays, Only : phir=>ubcp, nlphir=>nlbcp, onlphir=>onlbcp

  !Use nlphimainarrays, Only : phir=>phip, nlphir=>nlphip, onlphir&
  !     &=>onlphip
  USE  bcsetarrays, ONLY : gradphi=>ppr,fr,diffn 
  USE general_funcs
  USE dependent_functions
  USE interpolation
  USE initialize_flo
  IMPLICIT NONE 
  PRIVATE 
  
  REAL(prcn) ::  drm, rad, rad2, lc_time, lc_dist, da(2)
  REAL(prcn) :: sumforcepoint, sumforcenode,fphirmeanfluid(1)&
       &,fphirmeanloc(1),fphirmeanfluidloc(1)
  
  integer ::  is(ndim),iii(ndim),io(ndim), io2(ndim)
  integer, save :: unitflux, unitfor , unitfluxsum, unitsurftemp, unitnuss
  PUBLIC ::bcsetscal, interpolate_phidata
CONTAINS 
  SUBROUTINE bcsetscal(rks)!(phir,nlphir,onlphir)

    USE errormesgs
    USE general_funcs
    USE nlarrays, ONLY : fr1=>ur1
    IMPLICIT NONE 
    !REAL(prcn) , DIMENSION(:,:,:,:), INTENT(in) ::  phir, nlphir, onlphir

    !-----------------------------------------------------------------------
    !     local variables
    !
    INTEGER, INTENT(IN) :: rks
    INTEGER ::  ioffset,isp
    INTEGER :: i,j,k,l,m,n,sp, iphs
    REAL(prcn) ::  volume, area, upimod, utemp
    REAL(prcn) ::  tempt, areabdy(nbody), tempor,&
         & phisurface_tmp(nbody,nspmx), area_spec1, flux_tmp,&
         & nusselt_avg(4), source_mean(nspmx),&
         & tmp_nu_error_array(nerr_steps), sourcenuss(1) 


    !Diagnostic Variables 
    !REAL(prcn) ::  tempor1, tempor2

    CHARaCTER*80 :: FILENAME1, FILENAME2, FILENAME3, FILENAME4 
    CHARACTER(LEN=80) :: formfile
    LOGICAL:: filexist, isopen
#if PARALLEL
    REAL(prcn) :: frsumloc,frsum
    COMPLEX(prcn) :: fsumloc,fsum
#endif
!!$    unitfor = getnewunit(minunitno,maxunitno)
!!$
!!$    OPEN(unit=unitfor,file='forcing.diag',form='formatted'))

    formfile='formatted'
    IF(I_AM_NODE_ZERO.and.first_pass) THEN 
       FILENAME1 = TRIM(RUN_NAME)//'_scalflux'//'.dat'
       FILENAME2 = TRIM(RUN_NAME)//'_scalfluxsum'//'.dat'
       FILENAME3 = TRIM(RUN_NAME)//'_surftemp'//'.dat'
       FILENAME4 = TRIM(RUN_NAME)//'_nusselt_comps'//'.dat'
       CALL RUN_TIME_FILE_OPENER(unitflux,FILENAME1, formfile)
       CALL RUN_TIME_FILE_OPENER(unitfluxsum, FILENAME2, formfile)
       CALL RUN_TIME_FILE_OPENER(unitnuss, FILENAME4, formfile)
       IF(LUMPED_CAP) THEN 
          unitsurftemp = getnewunit(minunitno, maxunitno)
          CALL RUN_TIME_FILE_OPENER(unitsurftemp,FILENAME3,formfile)
       end IF

       !first_pass=.false.
    ENDIF



    drm = 1.d0 !its the value of the smearing length
    radm = 1.d0 !its the offset of the internal reversal point from the boundary point
    !-----------------------------------------------------------------------
    !     zero the force array

    Area_spec1 = (pi*(two*radbdy(1)*dx)**2)
    !phirmean(1:nspmx) = phirmean(1:nspmx) + fphirmean(1:nspmx)*dt
    !phigrad(1:nspmx) = zero
    DO sp=1,nspmx
       do k = 1, mz
          do j = 1, my
#if PARALLEL
             fr(0,j,k,sp)=zero 
#endif
             do i = 1, nx
                fr(i,j,k,sp) = zero
             ENDDO
#if PARALLEL
             fr(nx+1,j,k,sp)=zero 
#endif
          end do
       end do

    end DO
    fphirmean(:) = zero 
    fphirmeanfluid(:) = zero 
    fphirmeanloc(:) = zero
    fphirmeanfluidloc(:) = zero 

    CALL calc_diffn 

    sumforcepoint = 0.d0
    sumforcenode = 0.d0

    !    sum_flux_nm1(1:nbody) = zero

    flux_global(:) = zero 
    flux_global2(:) = zero 


    !-----------------------------------------------------------------------
    !     loop over all bodies

    DO m = 1, nbody !Loop over bodies
       iphs = 1!part_array(m)%iphs
       
       nbnd = phase_array(iphs)%nbnd
       nrpr = phase_array(iphs)%nrpr
       
       bndarray => phase_array(iphs)%bndpts
       
       da(1)=4.*pi*(radbdy(m)*dx)**2./float(nbnd)
       da(2)=4.*pi*(radibdy(m)*dx)**2./float(part_array(m)%nrpr_active)
       areabdy(m) = 4.*pi*(radbdy(m)*dx)**2
       !Write(*,*)'body = ', m 
       !     Calculate the surface area represented by each boundary point and 
       !     also do the same for each internal reversal point.

       !     Check if the forcing is conserved at the Boundary points
       flux_body(m,:) = zero
       flux_body2(m,:) = zero

       CALL calc_bndrev_data(m,rks)!phir,nlphir,onlphir, m)!force the boundary pts
       if(dorpr)  CALL  calc_inner_rev_data(m,rks)
       flux_global(:) = flux_global(:) + flux_body(m,:)*da(1)
#if !PARALLEL
       flux_global2(:) = flux_global2(:) + flux_body2(m,:)*da(1)
#endif

       !To get the forced data at the inner reversal points 

       !print*, '##############################################'
       !print*,'Check if forcing is conserved at the reversal pairs'
       !print*,'SPHERE NUMBER = ', m 
       !print*,'SUMFORCEPOINT AT rev pairs = ', sumforcepoint
       !print*,'SUMFORCENODE AT rev pairs = ', sumforcenode*dx**3.
       !Reset the surface temperature of the body based on lumped
       ! capacitance 

    END DO !end the body loop
    do isp = 1, nspmx
       GLOBAL_DOUBLE_SUM(fphirmeanloc(isp),fphirmean(isp),1,decomp_group)
    end do

    if(I_AM_NODE_ZERO)then
#if PARALLEL
       PRINT*,'flux_global = ', flux_global
#else
       PRINT*,'flux_global = ', flux_global, flux_global2
#endif
    end if
    !    if(I_AM_NODE_ZERO)PRINT*,'flux_global = ', flux_global!, flux_global2
    
    fphirmean(1:nspmx) = fphirmean(1:nspmx)/real((mx1*my*mz),prcn)
    fphirmeanfluid(1:nspmx) = fphirmeanfluid(1:nspmx)&
         &/real(count_fluid,prcn)

    if(I_AM_NODE_ZERO)WRITE(*,'(A25,2(2x,g17.8))')'FPHIRMEAN ', fphirmean!, fphirmeanfluid
    
    IF(setphimean) then 
       
       sourcesink(1:nspmx) = -fphirmean(1:nspmx)/(one-maxvolfrac)
       
       sourcenuss(:) = sourcesink(:) + fphirmeanfluid(:)
       
       if(I_AM_NODE_ZERO)WRITE(*,'(A30,2x,g17.8)')'NOT EVOLVING PHIRMEAN: SOURCE = ', sourcesink(1:nspmx)
       !WRITE(*,'(A30,2(2x,g17.8))')'fphirmean = ', fphirmean(1:nspmx), SUM(fr(1:mx1,:,:,1))/(mx1*my*mz)!sourcesink(1:nspmx)
       !phirmean does not evolve in this case
    ELSE

       
       sourcesink(:) = -cf*(phistream-phi_fluid_mean(:)) + gamma(:)*flux_global(:)/(voldom*(one&
            &-maxvolfrac))
       if(include_frmeanfluid) then 
          WRITE(*,'("INCLUDE_FRMEANFLUID IS TRUE")') 
          sourcesink(:) = sourcesink(:) - fphirmeanfluid(:)/(one - maxvolfrac)
       end if



       sourcenuss(:) = gamma(:)*flux_global(:)/(voldom*(one&
            &-maxvolfrac))

       phirmean(1:nspmx) = phirmean(1:nspmx) + (fphirmean(1:nspmx)&
            &+sourcesink(1:nspmx)*(one-maxvolfrac))*dt 

       if(I_AM_NODE_ZERO)then
          WRITE(*,'(A)')'-----------------EVOLVING PHIRMEAN----------------------'

          WRITE(*,'(A25,1x,g17.8)') 'UNSTEADY TERM:', cf*(phistream-phi_fluid_mean(1))*(one-maxvolfrac)*area_spec1&
               &/(6.d0*pi*maxvolfrac*gamma(1)*(phisurf-phistream))
          
          WRITE(*,'(A25,1x,g17.8)') 'DIFF TERM:', -(flux_global(1)*char_length*char_length)&
               &/(6.d0*voldom*maxvolfrac*(phisurf-phistream))
          
          WRITE(*,'(A25,1x,g17.8)') 'SOURCE TERM:', -(sourcesink(1)*(one-maxvolfrac)*area_spec1)&
               &/(6.d0*pi*maxvolfrac*gamma(1)*(phisurf-phistream))
          
          if(abs(frmeanfluid(1)).gt.zero) WRITE(*,'(A,2x,g17.8)') 'SOURCE/fphirmeanfluid', sourcesink(1)/(-fphirmeanfluid(:)/(one - maxvolfrac))
       end if
    end IF
    
    

    nusselt_avg(1) = -(sourcesink(1)*(one-maxvolfrac)*area_spec1)&
         &/(6.d0*pi*maxvolfrac*gamma(1)*(phisurf-phistream))
    nusselt_avg(2) = -(sourcenuss(1)*(one-maxvolfrac)*area_spec1)&
         &/(6.d0*pi*maxvolfrac*gamma(1)*(phisurf-phistream))
    nusselt_avg(3) = -(flux_global(1)*dchar*dchar)&
         &/(6.d0*voldom*maxvolfrac*(phisurf-phistream))

    !if(dorpr) then 
!!$    nusselt_avg(4) = -(flux_global2(1)*dchar*dchar)&
!!$         &/(6.d0*voldom*vol_frac1*(phisurf-phistream))
    !else
    !end if

    if(ABS(nusselt_avg(2)).gt.zero) then
       nu_error = ABS(nusselt_avg(1)-nu_old)/nusselt_avg(1)
    ELSE
       nu_error  = 1.d0
    end if


    if(rks.eq.itrmax) then 
       !Rearrange the ferror_array array so that the last entry is flushed out
       tmp_nu_error_array(1:nerr_steps) = nu_error_array(1:nerr_steps)
       nu_error_array(2:nerr_steps) = tmp_nu_error_array(1:nerr_steps-1)
       nu_error_array(1) = nu_error
       !PRINT*,'FERROR_A =', FERROR_ARRAY
       nu_error_hist = SUM(nu_error_array(1:nerr_steps))/nerr_steps
    end if

    if(I_AM_NODE_ZERO)WRITE(*,'(A25,6(2x,g14.7))')'NUs and ERROR(2)  = ', nusselt_avg(1:3), nu_error
    nu_old = nusselt_avg(1)
#if PARALLEL
    do isp = 1,nspmx
       VECSENDRECV(fr(nx+1,1,1,isp),1,twodrslice,toproc,1,diffn(0,1,1,isp),1,fromproc,1,decomp_group,status)
       VECSENDRECV(fr(0,1,1,isp),1,twodrslice,fromproc,1,diffn(nx+1,1,1,isp),1,toproc,1,decomp_group,status)
    end do
    
    DO isp = 1, nspmx 
       do j = 1, my 
          do k = 1, mz
             fr(1,j,k,isp) = fr(1,j,k,isp) + diffn(0,j,k,isp)
             fr(nx,j,k,isp) = fr(nx,j,k,isp) + diffn(nx+1,j,k,isp)
          end do
       end do
       !frsumloc = SUM(fr(1:nx,:,:,1))
       !GLOBAL_DOUBLE_SUM(frsumloc,frsum,1,decomp_group)
    end DO
    
    !if(I_AM_NODE_ZERO) WRITE(*,'(A25,3(2x,g17.8))')'fphirmean from domain=', frsum/(mx1&
    !     &*my*mz)
#endif

#if !PARALLEL    
    IF(xperiodic)  fr(mxf, :, :,1:nspmx) = fr(1,:,:,1:nspmx)
#endif
    !PRINT*,'MAXVOLFRAC = ', maxvolfrac
    !source_mean = zero 
    DO i=1,nx !mx1             !loop over planes
       DO sp =1, nspmx   
          do j = 1, my 
             do k = 1, mz

                fr(i,j,k,sp) = fr(i,j,k,sp) - fphirmean(sp)
                !sumforcenode=sumforcenode+fr(i,j,k,sp)
                if(fluid_atijk(i,j,k)) then 

                   fr(i,j,k,sp)  = fr(i,j,k,sp)+sourcesink(sp)*maxvolfrac
                else
                   fr(i,j,k,sp) = fr(i,j,k,sp) - sourcesink(sp)*(one&
                        &-maxvolfrac)
	!	   source_mean(sp) = source_mean(sp)+fr(i,j,k,sp)!-sourcesink(sp)*(one-maxvolfrac)
                end if
                
             end do
          end do
          CALL ff2rc(fr(i,:,:,sp),ffphi(i,:,:,sp))
          !ffphi(i,1,1,sp) = czero
       ENDDO
       !       Print*,'force  = ', ff(i,1,1,1),i 
    ENDDO
#if !PARALLEL
    ffphi(mx,:,:,1:nspmx) = ffphi(1,:,:,1:nspmx)
#endif

#if PARALLEL
    !fsumloc = SUM(ff(1:nx,1,1,1))
    !GLOBAL_COMPLEX_SUM(fsumloc,fsum,1,decomp_group)
    !if(I_AM_NODE_ZERO)WRITE(*,'(A25,10(2x,g17.8))')'FLO: AVG FLUC FORCE  = ',fsum/(mx1)
#else
    WRITE(*,'(A25,10(2x,g17.8))')'AVG FLUC FORCING = '&
         &,SUM(ffphi(1:mx1,1,1,1))/(mx1)
#endif


    if (I_AM_NODE_ZERO)then
       IF(rks.eq.itrmax) then 

          WRITE(unitflux,'(400(2x,e20.12))') (t)/t_conv, t/t_vis, t/t_diff,&
               & (-flux_body(i,1)/(phisurf-phistream), i=1,nbody) , -flux_global(1)/(phisurf-phistream)
          
          WRITE(unitnuss,'(20(2x,e20.12))')  nusselt_avg(3),nusselt_avg(1), ABS(nusselt_avg(3)-nusselt_avg(1))/(SMALL_NUMBER + nusselt_avg(3))
          
          WRITE(unitfluxsum,'(20(2x, e20.12))') t/t_conv, t/t_vis, t&
               &/t_diff, nusselt_avg(3),nusselt_avg(1), ABS(nusselt_avg(3)-nusselt_avg(1))/nusselt_avg(3), nu_error, phi_fluid_mean(1)/phistream
       end IF

    end if
  END SUBROUTINE bcsetscal

  SUBROUTINE calc_bndrev_data(m,rks)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: rks, m

    REAL(prcn) ::  dfll(nspmx)

    REAL(prcn) ::  xl(ndim), phil(nspmx),philo(nspmx),philo2(nspmx),phili(nspmx)
    REAL(prcn) ::  nlphil(nspmx),onlphil(nspmx),xltemp
    INTEGER :: sp, i, j, k, n, l 
    REAL(prcn) :: rad, tempor(nspmx), gradphibnd(3,nspmx), normal(3),&
         & sourcell(nspmx),fluxloc(nspmx)

    INTEGER :: ib, ie, jb, je, kb, ke, onew, ii, jj, kk, phicelltemp

    Integer :: pcell(3)
    LOGICAL :: phiterm

    frombnd = .TRUE.
    fromrpr = .FALSE.
    fluxloc = zero

    DO  10 l=1,nbnd
       
       rad = zero
       phil(:)=zero
       nlphil(:)=zero
       onlphil(:)=zero
       dfll(:) = zero
       sourcell(:) = zero 

       DO n=1,ndim
          
          xl(n)=xc(m,n)+ bndarray(n,l)*radbdy(m)
          
          
          rad=rad+(bndarray(n,l)*radbdy(m))**2.0
          if(xl(n).lt.zero) then 
             is(n) = int(xl(n)-1)
          else 
             is(n) = int(xl(n))
          end if

       ENDDO

       rad = SQRT(rad)
       normal(1:3) = bndarray(1:3,l)

#if PARALLEL
       
       xltemp = xl(1)
       phicelltemp = is(1)
       if(.not.CELL_IN_PROC(phicelltemp))then
          WEST_PERIODIC_IMAGE(is(1),phicelltemp,xl(1),xltemp)
          EAST_PERIODIC_IMAGE_MOD(is(1),phicelltemp,xl(1),xltemp)
          if(.not.CELL_IN_PROC(phicelltemp))goto 2600
       end if
       phiterm = CELL_IN_VEL_GRID(phicelltemp)
       
       xl(1) = xltemp
       is(1) = phicelltemp
#else
       phiterm = .TRUE.
#endif

       call interpolate_phidata(is,xl, ib&
            &,ie,jb,je,kb,ke,phil,nlphil,onlphil,dfll, 1,onew) 
            
       DO sp = 1,nspmx 
          !surf_scal_value(m,l,sp) = phil(sp)
          tempor(sp) = zero
          if(phiterm)tempor(sp) = tempor(sp) +cf*(phil(sp)-phisurfall(m,sp))
          tempor(sp) = tempor(sp)- coef(rks,3)*nlphil(sp)-coef(rks,4)*onlphil(sp) 
          tempor(sp) = tempor(sp) -(coef(rks,1)+coef(rks,2))*dfll(sp)
          tempor(sp) = tempor(sp)*da(1)*drm*dx
          sumforcepoint = sumforcepoint + tempor(sp)/dx**3.d0
       ENDDO
       
       do k = 1, onew
          kk = kb+k-1
          if(kk.lt.1) kk = mz+kk
          if(kk.gt.mz) kk = kk-mz 
          
          do j = 1, onew
             jj = jb+j-1
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             
             do i = 1, onew
                DO sp = 1, nspmx
                   ii = ib+i-1
#if !PARALLEL
                   if(ii.lt.1) ii = mxf+ii-1
                   if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
                   LOCAL_INDEX(ii)
#if PARALLEL
                   if(ii.lt.0)then
                      PRINT*,'BND :ii is less than 0', myid,ii,xl(1)
                   else if(ii.gt.nx+1) then
                      PRINT*,'BND :ii is gt than nx+1', myid,ii,xl(1)
                   end if
#endif
                   IF(DOBND) then 
                      fr(ii,jj,kk,sp)=fr(ii,jj,kk,sp)+(weightp(i,j,k)&
                           &*tempor(sp))/dx**3.  
                      fphirmeanloc(sp) = fphirmeanloc(sp) + weightp(i,j,k)*tempor(sp)/dx**3. 
                      if(fluid_atijk(ii,jj,kk)) fphirmeanfluidloc(sp) =&
                           & fphirmeanfluidloc(sp) + (weightp(i,j,k)&
                           &*tempor(sp))/dx**3.   
                   END IF

                   gradphisten(i,j,k,:,sp) = gradphi(ii,jj,kk,:)

                   !if(fluid_atijk(ii,jj,kk)) fphirmeanfluidloc(sp) =&
                   !    & fphirmeanfluidloc(sp) + (weightp(i,j,k)&
                   !   &*tempor(sp))/dx**3.   
                   
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    
    do sp = 1,nspmx
       
       do n = 1, ndim 
          gradphibnd(n,sp) = &
               & array_dot_product(gradphisten(1:onew,1:onew&
               &,1:onew,n,sp),weightp(1:onew,1:onew,1:onew))
       end do
       fluxloc(sp) = fluxloc(sp)+&
            & array_dot_product(gradphibnd(1:ndim,sp)&
            &,normal(1:ndim)) 
       
!!$          flux_body(m,sp) = flux_body(m,sp)+&
!!$               & array_dot_product(gradphibnd(1:ndim,sp)&
!!$               &,normal(1:ndim)) 
    end do

       !RINT*,'flux = ', flux_body(m,1)
          
2600 continue
      
10  ENDDO    !close loop over all boundary points
    do sp=1,nspmx
       GLOBAL_DOUBLE_SUM(fluxloc(sp),flux_body(m,sp),1,decomp_group)
    end do

   frombnd = .FALSE.

  END SUBROUTINE calc_bndrev_data

  SUBROUTINE calc_inner_rev_data(m,rks)

    IMPLICIT NONE 

    ! REAL(prcn) , DIMENSION(:,:,:,:), INTENT(in) ::  phir, nlphir, onlphir
    REAL(prcn) ::  dfll(nspmx),  xlo(ndim),xli(ndim), xlo2(ndim)
    INTEGER, INTENT(IN) :: rks

    REAL(prcn) ::  xl(ndim), phil(nspmx),force_dist(3), tmppa, force_fl(3)
    REAL(prcn) ::  philo(nspmx),philo2(nspmx),phili(nspmx)
    REAL(prcn) ::  nlphil(nspmx),onlphil(nspmx), xltemp
    INTEGER, INTENT(in) :: m
    Integer :: pcell(3)
    INTEGER :: sp, i, j, k,l, count_fl, count_so, n
    REAL(prcn) :: rad, tempor(nspmx), sourcell(nspmx)
    INTEGER :: ib, ie, jb, je, kb, ke, onew, ii, jj, kk,rprcountloc,&
         & rprcount,rprcom,rprcomloc, phicelltemp
    LOGICAL :: phiterm
#if PARALLEL
    INTEGER :: FOCUS_POINT, FOCUS_PARTICLE
    FOCUS_POINT = -1
    FOCUS_PARTICLE = -1
#endif

    frombnd = .FALSE.
    fromrpr = .TRUE.

   
    rprcount = 0
    rprcountloc = 0
    rprcom = 0
    rprcomloc = 0
    DO l=1,nrpr
       if(.NOT.PART_ARRAY(m)%if_rev(L)) GOTO 666
       rad = zero
       rad2 = zero

       DO sp=1,nspmx
          phili(sp)=zero
          philo(sp)=zero
          nlphil(sp)=zero
          onlphil(sp)=zero
          dfll(sp) = zero
          philo2(sp)=zero
          sourcell(sp) = zero 
       ENDDO

       DO 20  n=1,ndim

          !     location of internal points

          xli(n)=xc(m,n)+ bndarray(n,l)*radibdy(m)

          !location of external points

          xlo(n)=xc(m,n)+ bndarray(n,l)*radobdy(m)

          xlo2(n) = xc(m,n)+ bndarray(n,l)*rado2bdy(m)

          rad=rad+(bndarray(n,l)*radobdy(m))**2.
          
          rad2=rad2+(bndarray(n,l)*rado2bdy(m))**2.
#if !PARALLEL
          if(xlo2(n).lt.zero) then 
             io2(n) = int(xlo2(n)-1)
          else 
             io2(n) = int(xlo2(n))
          end if
#endif

          if(xlo(n).lt.zero) then 
             io(n) = int(xlo(n)-1)
          else 
             io(n) = int(xlo(n))
          end if
          
          if(xli(n).lt.zero) then 
             iii(n) = int(xli(n)-1)
          else 
             iii(n) = int(xli(n))
          end if
       
20     ENDDO
       
       rad=dsqrt(rad)
       rad2=dsqrt(rad2)

#if PARALLEL
       phicelltemp = iii(1)
       xltemp  = xli(1)
       if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
          PRINT*,' INNER REV PT = ', xli(1),phicelltemp, myid, m
          PRINT*,' EXTERNAL REV PT = ', xlo(1), myid,m
       end if
       
       if(.not.CELL_IN_PROC(phicelltemp))then
          WEST_PERIODIC_IMAGE(iii(1),phicelltemp,xli(1),xltemp)
          EAST_PERIODIC_IMAGE_MOD(iii(1),phicelltemp, xli(1), xltemp)
          if(.not.CELL_IN_PROC(phicelltemp))then
             if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
                PRINT*,' INNER REVERSAL PT = ', xli(1),phicelltemp, myid, m
                
                !PARALLEL_FINISH()
                !STOP
             end if
             goto 666
          end if
       end if
       if(EAST_NO_MANS_LAND(iii(1)).or.EAST_NO_MANS_LAND(phicelltemp)) then 
          phiterm = .not.CONCAVE(xli,1,m)
       else if(WEST_NO_MANS_LAND(iii(1)).or.WEST_NO_MANS_LAND(phicelltemp))then
          phiterm = CONCAVE(xli,1,m)
       else
          phiterm = .TRUE.
       end if
       iii(1) = phicelltemp
       xli(1) = xltemp

       if(phiterm)then
          phicelltemp = io(1)
          xltemp = xlo(1)
          if(.not.RPR_CELL_IN_PROC(phicelltemp))then
             WEST_PERIODIC_IMAGE(io(1),phicelltemp, xlo(1),xltemp)
             EAST_PERIODIC_IMAGE_MOD(io(1),phicelltemp,xlo(1),xltemp)
             if(.not.RPR_CELL_IN_PROC(phicelltemp))then
                if(I_AM_NODE_ZERO)then
                   if(phicelltemp.eq.mxf-3)then
                      phicelltemp = phicelltemp-(mxf-1)+1
                      xltemp = xltemp-(mxf-1)
                   endif
                else
                   PRINT*,' ERROR WITH EXTERNAL POINT IN THIS PROCESSOR : ', myid, m, l, xlo(1), phicelltemp,io(1),xli(1)
                end if
             end if
          end if
          io(1) = phicelltemp
          xlo(1) = xltemp
       end if

       if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
          PRINT*,' PHITERM = ', l,phiterm, myid, m
          PARALLEL_FINISH()
          STOP
       endif
#else
       phiterm = .TRUE.
#endif

       if(phiterm)then
          call interpolate_phidata(io,xlo, ib&
               &,ie,jb,je,kb,ke,philo,nlphil,onlphil,dfll, 0,onew) 
       end if

#if !PARALLEL
       call interpolate_phidata(io2,xlo2, ib&
            &,ie,jb,je,kb,ke,philo2,nlphil,onlphil,dfll, 0,onew)
#endif

       call interpolate_phidata(iii,xli,ib&
            &,ie,jb,je,kb,ke,phili,nlphil,onlphil,dfll, 1,onew) 

#if !PARALLEL
       DO sp = 1, nspmx
          !Nu3(m,l,sp)=-(phisurfall(m,sp) - philo(sp))/  &
          !    & ((-radbdy(m)+radobdy(m))*dx)

          flux_body2(m,sp)= flux_body2(m,sp) + (-three*phisurfall(m,sp)&
               &+four*philo(sp)-philo2(sp))/(two*(-radbdy(m)&
               &+radobdy(m))*dx) 
          !    Print*,'Nu=', Nu3(m,l,sp), maxval(-Nu3(m,:,sp))/minval(-Nu3(m,:,sp))
       END DO
#endif
       if(phiterm)then
          DO sp=1,nspmx
             
             philo(sp) = phisurfall(m,sp)*(radobdy(m)-radibdy(m))&
                  &-philo(sp)*(radbdy(m)-radibdy(m))
             philo(sp) = philo(sp)/(radobdy(m)-radbdy(m))
             
             phil(sp)=cf*(-philo(sp)+phili(sp))

          ENDDO
       end if

       !----------------------------------------------------------
       !     Interpolate the forcing using the interpolating
       ! factor on the boundary.

       !NOTE: IN ORDER TO PUT CORRECT FORCING ON THE GRID THE LAST
       ! INTERPOLATION SHUD HAVE BEEN PERFORMED ON THE INNER POINT OR

       DO sp = 1, nspmx
          tempor(sp) = zero
          if(phiterm) tempor(sp) = tempor(sp) + phil(sp)
          tempor(sp) = tempor(sp) - coef(rks,3)*nlphil(sp) -coef(rks&
               &,4)*onlphil(sp)
          tempor(sp) = tempor(sp) - (coef(rks,2)+coef(rks,1))*dfll(sp)
          tempor(sp) = tempor(sp)*da(2)*drm*dx  

          !fphirmean(sp) = fphirmean(sp) + tempor(sp)/dx**3.d0
          
          sumforcepoint = sumforcepoint + tempor(sp)/dx**3.d0
       ENDDO
       count_fl = 0 
       force_fl = zero 

       do k = 1, onew
          kk = kb+k-1
          if(kk.lt.1) kk = mz+kk
          if(kk.gt.mz) kk = kk-mz 
          
          do j = 1, onew
             jj = jb+j-1
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             
             do i = 1, onew
                
                ii = ib+i-1
#if !PARALLEL
                if(ii.lt.1) ii = mxf+ii-1
                if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
                LOCAL_INDEX(ii)
                
                if(ii.lt.0)then
                   PRINT*,'ii is less than 0', myid,ii,xli(1)
                else if(ii.gt.nx+1) then
                   PRINT*,'ii is gt than nx+1', myid,ii,xli(1)
                end if
                if(fluid_atijk(ii,jj,kk)) then 
                   count_fl = count_fl+1
                   DO sp = 1, nspmx
                      tmppa = weightp(i,j,k)*tempor(sp)/dx**3.  
                      force_fl(sp) = force_fl(sp)+tmppa
                   end DO
                end if
             end do
          end do
       end do
       count_so = onew*onew*onew - count_fl
       force_dist(1:nspmx) = force_fl(1:nspmx)/real(count_so,prcn)
       do k = 1, onew
          kk = kb+k-1
          if(kk.lt.1) kk = mz+kk
          if(kk.gt.mz) kk = kk-mz 
          
          do j = 1, onew
             jj = jb+j-1
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             
             do i = 1, onew
                ii = ib+i-1
#if !PARALLEL
                if(ii.lt.1) ii = mxf+ii-1
                if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
                
                LOCAL_INDEX(ii)
                
                if(include_frmeanfluid)  then 
                   if(fluid_atijk(ii,jj,kk)) then
                      DO sp = 1, nspmx
                         fphirmeanfluidloc(sp) =&
                              & fphirmeanfluidloc(sp)  +weightp(i,j,k)&
                              &*tempor(sp)/dx**3.
                      End DO
                   endif
                   DO sp = 1, nspmx
                      fr(ii,jj,kk,sp)=fr(ii,jj,kk,sp)+weightp(i,j,k)&
                           &*tempor(sp)/dx**3.
                   End DO
                else
                   if(.not.fluid_atijk(ii,jj,kk)) then
                      DO sp = 1, nspmx
                         fr(ii,jj,kk,sp)=fr(ii,jj,kk,sp)+weightp(i,j,k)&
                              &*tempor(sp)/dx**3.  + force_dist(sp)
                      End DO
                   end if
                end if
                
                fphirmeanloc(1:nspmx) = fphirmeanloc(1:nspmx) + weightp(i,j,k)*tempor(1:nspmx)/dx**3. 
                
                !if(fluid_atijk(ii,jj,kk)) fphirmeanfluidloc(sp) =&
                     !    & fphirmeanfluidloc(sp) + (weightp(i,j,k)&
                !   &*tempor(sp))/dx**3.   
             ENDDO
          ENDDO
       ENDDO


666    continue
    ENDDO ! loop over reversal points
    fromrpr = .FALSE.


  END SUBROUTINE calc_inner_rev_data


  SUBROUTINE calc_diffn 
    USE nlarrays, ONLY : ff1=>uf1, ff2=>uf2, ff3=>uf3, ff4=>uf11,fr1&
         &=>ur1,fr2=>ur2,fr3=>ur3 
    IMPLICIT NONE 
    INTEGER :: sp, i, j, k, ioffset
    Real(prcn) :: sumdiff1, sumdiff2, sumdiff3
    sumdiff1 = zero
    sumdiff2 = zero
    sumdiff3 = zero


    DO sp =1, nspmx ! Initialize the loop for all the species (nspmx)
#if PARALLEL
       diffn(0,:,:,sp) = zero
       gradphi(0,:,:,:)= zero
#endif
       DO i=1,nx !mxf-1 !Was 2,mxf previously for reasons unknown RG 02/28/06
          ioffset=i
          DO k=1,mz
             DO j=1,my2
                if(xperiodic) then
#if PARALLEL
                   ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                        &,k,sp)+phif(ioffset+1,j,k,sp)) 
#else                   
                   if(ioffset.eq.1) THEN
                      ff1(j,k)=(1./dx2)*(phif(nx,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+ phif(ioffset+1,j,k,sp)) 
                   else if(ioffset.eq.nx) THEN
                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(1,j,k,sp))

                   else

                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(ioffset+1,j,k,sp)) 
                   endif
#endif
	        else
                   ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.&
                        &*phif(ioffset,j,k,sp)+phif(ioffset+1,j,k,sp)) 
                endif
                ff1(j,k)=ff1(j,k)-w2(j,k)*phif(ioffset,j,k,sp)
                
                ff1(j,k)=ff1(j,k)*gamma(sp)
                
                !end of diffusion terms calculation 
                !beginnig of grad phi tems calculation
#if PARALLEL
                ff2(j,k) = (phif(i+1,j,k,sp)-phif(i-1,j,k,sp))/(two&
                     &*dx)
#else
		
                   if(i.eq.1) then
                      ff2(j,k) = (phif(2,j,k,sp)-phif(mxf-1,j,k,sp))/(two*dx) 

                   else if(i.eq.(nx))then
                      ff2(j,k) = (phif(1,j,k,sp)-phif(i-1,j,k,sp))/(two*dx) 
                   ELSE
                      ff2(j,k) = (phif(i+1,j,k,sp)-phif(i-1,j,k,sp))/(two*dx)
                   end if
#endif
                   ff3(j,k)=phif(i,j,k,sp)*wy(j)  !starts at foffset+1
                   ff4(j,k)=phif(i,j,k,sp)*wz(k)
                   

                ENDDO
             ENDDO

          CALL ff2cr(ff1(:,:),diffn(i,:,:,sp))
          
          CALL ff2cr(ff2(:,:),gradphi(i,:,:,1))
          CALL ff2cr(ff3(:,:),gradphi(i,:,:,2))
          CALL ff2cr(ff4(:,:),gradphi(i,:,:,3))
       
       END DO
#if PARALLEL
       diffn(nx+1,:,:,sp) = zero
       gradphi(nx+1,:,:,:)= zero
#endif
    ENDDO
    !write(90,'(12(2x,e20.10))') t, sumdiff1/(mxf*my*mz), sumdiff2/(mxf&
    !     &*my*mz), sumdiff3/(mxf*my*mz)
  END SUBROUTINE calc_diffn



  subroutine interpolate_phidata(pc,pos,ib,ie,jb,je,kb,ke,ul&
       &,nll,onll,dfll,flag, onew)

    USE general_funcs
    USE interpolation

    implicit none 
    Integer, intent(in) :: pc(3)
    !REAL(prcn), DIMENSION(:,:,:,:), INTENT(in) ::  ur, nlr, onlr
    Integer, Intent(in) :: flag
    Real(prcn), Dimension(:), Intent(in) :: pos
    Integer, INtent(out) :: ib, ie, jb,je,kb,ke, onew
    Real(prcn), Intent(out), Dimension(:) :: ul,nll,onll,dfll
    Integer :: i, j,k, ii,jj,kk, n

    CALL set_interpolation_stencil(pc,ib,ie,jb,je,kb,ke&
         &,interp_scheme, onew) 
    if((ib.lt.1.or.ie.gt.mxf).AND.(.not.xperiodic)) Print*,'Error in i ....',ib,ie,pc,pos

    do k = 1, onew
       do j = 1, onew
          do i = 1, onew
             ii = ib+i-1
             jj = jb+j-1
             kk = kb+k-1
             gstencil(i,j,k,1) = ib+(i-1)
             gstencil(i,j,k,2) = jb+(j-1)
             gstencil(i,j,k,3) = kb+(k-1)
#if !PARALLEL
             if(ii.lt.1) ii = mxf+ii-1
             if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             if(kk.lt.1) kk = mz+kk
             if(kk.gt.mz) kk = kk-mz
!!$   
             LOCAL_INDEX(ii)
#if PARALLEL
             if(flag.eq.1)then
                if(ii.lt.0) then 
                   PRINT*,'phidata :ERROR WEST SIDE', myid,pos(1), frombnd&
                        &,fromrpr
                   PARALLEL_FINISH()
                   STOP
                else if(ii.gt.nx+1)then
                   PRINT*,'phidata :ERROR EAST SIDE', myid,pos(1),&
                        & frombnd,fromrpr
                   PARALLEL_FINISH()
                   STOP
                end if
             else
                if(ii.lt.-1) then 
                   PRINT*,'phidata external :ERROR WEST SIDE', myid,pos(1), frombnd&
                        &,fromrpr
                   PARALLEL_FINISH()
                   STOP
                else if(ii.gt.nx+2)then
                   PRINT*,'phidata external :ERROR EAST SIDE', myid,pos(1),&
                        & frombnd,fromrpr
                   PARALLEL_FINISH()
                   STOP
                end if

             end if
#endif
             vsten(i,j,k,1:nspmx) = phir(ii,jj,kk,1:nspmx)
             if(flag.eq.1) then 
                nlsten(i,j,k,1:nspmx) = nlphir(ii,jj,kk,1:nspmx)
                onlsten(i,j,k,1:nspmx) = onlphir(ii,jj,kk,1:nspmx)
                dfsten(i,j,k,1:nspmx) = diffn(ii,jj,kk,1:nspmx)
             end if

          end do
       end do
    end do
    CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:3),&
         & vsten(1:onew,1:onew,1:onew,1:nspmx),pos(1:ndim),ul(1:nspmx),onew,&
         & interp_scheme,weightp) 
    if(flag.eq.1) then 
       do n = 1, nspmx 
          nll(n) =  array_dot_product(nlsten(1:onew,1:onew,1:onew,n&
               &),weightp(1:onew,1:onew,1:onew)) 
          onll(n)=  array_dot_product(onlsten(1:onew,1:onew,1:onew,n&
               &),weightp(1:onew,1:onew,1:onew))
          dfll(n)=  array_dot_product(dfsten(1:onew,1:onew,1:onew,n&
               &),weightp(1:onew,1:onew,1:onew))
       end do

    end if
    !Print*, 'nll in interpdata =', nll(1), nlsten
  end subroutine interpolate_phidata

END MODULE bc_scalar
  
