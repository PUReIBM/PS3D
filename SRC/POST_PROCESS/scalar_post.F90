MODULE scalar_post
  USE precision 
  USE constants 
  USE scalar_data
  Use nlmainarrays, Only : phir=>ubcp, nlphir=>nlbcp
  USE general_funcs
  USE postproc_funcs
  USE post_global_data
  USE fftw_interface
  !
  USE  bcsetarrays, ONLY : gradphi=>ppr, diffn, fr

  implicit none 
  Private 
  REAL(prcn) :: axial_direction(3), sum_nlterm
  REAL(prcn), dimension(:,:,:), Allocatable :: flux_body_components
  REAL(prcn), dimension(:), Allocatable :: flux_global_parallel,&
       & flux_global_perp
  REAL(prcn) :: sourcesink2
  Public :: scalar_pp
CONTAINS
  subroutine  scalar_pp(imis)
    USE nlphi 
    Use nlarrays , only : ur1,uf1
    implicit none 
    INTEGER, Intent(in) :: imis
    LOGICAL, SAVE :: routine_called = .false. 
    Integer :: i, j, k, count_tmp, runitno, nrbins, ibin, isp
    INTEGER, SAVE :: count_routine
    CHARACTER*80 :: SCALANDU_FOR, THETAVSNU3,SCALANDU,SCALONLY, junk_char, filename1, stat 
    LOGICAL :: filexist , isopen
    real(prcn), dimension(:), allocatable :: phi_corr, rad_bin, conf_l95, conf_r95

    REAL(PRCN) :: total_flux, deficit
    IF(routine_called) then 
       DEALLOCATE(flux_body_components, flux_global_parallel, flux_global_perp)
    end IF
    
    WRITE(*,*) 'nspmx in scalar_pp = ', nspmx
    routine_called = .true.
    nrbins = NINT(real(my,prcn)/2.d0) + 1

    axial_direction(1:3) = uchar(1:3)/ucharmod

    WRITE(*,'(A,3(2x,g17.8))') 'axial direction = ', axial_direction(1:3)

    ALLOCATE(phi_corr(nrbins), rad_bin(nrbins), conf_l95(nrbins), conf_r95(nrbins))

    CALL nlphi_graduphi(1)

    CALL calc_diffn
    call calc_flux(imis) 

    nlcontrib(imis, :) = zero 
    diffcontrib(imis, :) = zero 
    fphicontrib(imis, :) = zero 
    fphicontrib_solid(imis, :) = zero 
    
    DO k = 1, mz
       do j = 1, my
          do i = 1, mx1
             do isp = 1, nspmx
                IF(fluid_atijk(i,j,k)) then 
                   nlcontrib(imis, isp) = nlcontrib(imis, isp) + nlphir(i,j,k,isp)

                   diffcontrib(imis, isp) = diffcontrib(imis, isp) + SUM(diffn(i,j,k,1:3))
                   fphicontrib(imis,isp) = fphicontrib(imis, isp) +&
                        & fr(i,j,k,isp)
                ELSE
                   fphicontrib_solid(imis,isp) = fphicontrib_solid(imis, isp) +&
                        & fr(i,j,k,isp)
                end IF
                
             end do
          end do
       end do
    end DO
    total_flux = flux_global_parallel(1) + flux_global_perp(1)
    

    WRITE(*,'(A, 2x, g17.8)')'TOTAL FLUX = ', total_flux
    sourcesink2 =  gamma(1)*total_flux/(voldom*(one&
         &-maxvolfrac))

    CALL compute_budget_along_flow(imis)

    deficit = -cf*(phistream-phi_fluid_mean(1))
    WRITE(*,'(A,2(2x,g29.15))') 'phis - phif = ',(phistream-phi_fluid_mean(1)), -cf
    WRITE(*,'(A, 2x, g17.8)')'SOURCE FROM ONLY FLUX = ', sourcesink2
    WRITE(*,'(A,2( 2x, g17.8))')'SOURCE - source2 = ', sourcesink(1)-sourcesink2, deficit

    nlcontrib(imis, :) = nlcontrib(imis, :)/real(count_fluid, prcn)
    diffcontrib(imis, :) = diffcontrib(imis,:)/real(count_fluid,&
         & prcn) 
    fphicontrib(imis, :) = fphicontrib(imis,:)/real(count_fluid, prcn)
    fphicontrib_solid(imis, :) = fphicontrib_solid(imis,:)/real(count_solid, prcn)

    WRITE(*,'(A30,2x,g17.8)') 'nl = ', -nlcontrib(imis,1)
    WRITE(*,'(A30,2x,g17.8)') 'nlterm at bnd = ', -sum_nlterm/voldom
    WRITE(*,'(A30,2x,g17.8)') '-diff = ', -diffcontrib(imis,1)
    WRITE(*,'(A30,2x,g17.8)') '-fr = ', -fphicontrib(imis,1)
    WRITE(*,'(A30,2x,g17.8)') 'nl-diff-fr = ', -nlcontrib(imis,1) - diffcontrib(imis,1) - fphicontrib(imis,1)
    WRITE(*,'(A30,2(2x,g17.8))') 'fr solid, fr fluid/source = ', &
         &-fphicontrib_solid(imis,1), -fphicontrib(imis,1)/sourcesink(1)


    WRITE(*,'(A30,2(2x,g17.8))') 'sources = ', sourcesink(1), sourcesink2
    !CALL scalar_two_point_correlation_pencil(phir(1:mx,1:my,1:mz,1:nspmx),phi_corr,rad_bin,conf_l95,conf_r95, nrbins)

!!$    runitno = getnewunit(minunitno, maxunitno)
!!$
!!$    open(unit=runitno, file = TRIM(RUN_NAME)//'_PHI_TWOPT_CORR.dat', form="formatted", status = "unknown")
!!$    do ibin=1,nrbins
!!$     write(runitno,'(10(2x,f12.8))') rad_bin(ibin), rad_bin(ibin)/doml(1), phi_corr(ibin), conf_r95(ibin)-phi_corr(ibin)
!!$    end do
!!$    
!!$    
!!$    close(runitno,status='keep')




    DEALLOCATE(phi_corr,rad_bin, conf_l95, conf_r95)




  end subroutine scalar_pp

  SUBROUTINE   compute_budget_along_flow(imis)
    Implicit none 
    INTEGER, Intent(in) ::imis
    REAL(prcn) :: nlpldiff, nlonly, diffonly, diff1, diff2, diff3, norm_quant, area_face, fluid_area_frac , phirpl, norm_quant1, fphir
    integer :: count_fluid_plane, i, j, k

    area_face = real(my,prcn)*real(mz,prcn)
    norm_quant1 = area_face
    
    
    CALL screen_separator(100, 'B')
    XLOOP: DO i = 1,mx1
       nlpldiff = zero
       count_fluid_plane = 0
       nlonly = zero
       diff1 = zero
       diff2 = zero
       diff3 = zero
       phirpl = zero
       fphir  = zero
       !REMEMBER NON-LINEAR TERM IS MULTILPIED BY MINUS SIGN IN NLPHI. THEREFORE, HERE NL IS AGAIN MULTIPLIED BY MINUS SIGN. NOT A BUG AND DONT GO RUNNING TO YOUR MAJOR PROFESSOR THAT I SCREWED UP ;). Rahul Garg 12/10/08
       DO j = 1, my
          DO k = 1, mz
             IF(fluid_atijk(i,j,k)) then 
                nlpldiff = nlpldiff - SUM(diffn(i,j,k,1:3)) - nlphir(i,j,k,1)-fr(i,j,k,1)
                nlonly = nlonly + nlphir(i,j,k,1)
                diff1 = diff1 + (diffn(i,j,k,1))
                diff2 = diff2 + (diffn(i,j,k,2))
                diff3 = diff3 + (diffn(i,j,k,3))
                phirpl = phirpl + phir(i,j,k,1)
                fphir = fphir + fr(i,j,k,1) 
                count_fluid_plane = count_fluid_plane + 1
             end IF
          end DO
       end DO
!!$       scal_frac_nl(i,1,imis) = nlonly/nlpldiff
!!$       scal_frac_diff(i,1,imis) = diff1/nlpldiff
!!$       scal_frac_diff(i,2,imis) = diff2/nlpldiff
!!$       scal_frac_diff(i,3,imis) = diff3/nlpldiff
       solid_area_frac(i, imis) = one - REAL(count_fluid_plane,prcn)/REAL(my*mz,prcn)
       fluid_area_frac = REAL(count_fluid_plane,prcn)/REAL(my*mz,prcn)

       
       norm_quant = norm_quant1*sourcesink(1)*fluid_area_frac

       scal_frac_nl_diff(i,1,imis) = nlpldiff/norm_quant

       scal_frac_nl(i,1,imis) = -nlonly/norm_quant
       scal_frac_diff(i,1,imis) = -diff1/norm_quant
       scal_frac_diff(i,2,imis) = -(diff2+diff3)/norm_quant
       scal_frac_diff(i,3,imis) = -diff3/norm_quant
       scal_frac_fphi(i,1,imis) = -fphir/norm_quant
       write(*,'(7(2x,g17.8))')  REAL(count_fluid_plane,prcn)/REAL(my*mz,prcn), scal_frac_nl(i,1,imis), scal_frac_diff(i,1,imis), scal_frac_diff(i,2,imis), scal_frac_fphi(i,1,imis), scal_frac_nl_diff(i,1,imis)!,sourcesink(1)*fluid_area_frac
    end DO XLOOP


!!$    WRITE(*,'("FROM BUDGET:","nl = ", g17.8,/,"-diff = ",g17.8,/,"-fr &
!!$         &= ", g17.8,/, "nl-d&
!!$         &iff-fr = ",g17.8)') SUM(scal_frac_nl(1:mx1,1,imis))&
!!$         &*area_face/real(count_fluid,prcn),&
!!$         & (SUM(scal_frac_diff(1:mx1,1,imis))&
!!$         &+SUM(scal_frac_diff(1:mx1,2,imis)))*area_face&
!!$         &/real(count_fluid,prcn), &
!!$         & SUM(scal_frac_fphi(1:mx1,1,imis))*area_face&
!!$         &/real(count_fluid,prcn) ,SUM(scal_frac_nl_diff(1:mx1,1&
!!$         &,imis))*area_face/real(count_fluid,prcn)

    
    CALL screen_separator(100, 'B')
  end SUBROUTINE compute_budget_along_flow

  subroutine calc_flux(imis)
    USE dependent_functions
    USE bc_scalar 
    implicit none
    Integer, Intent(in) :: imis
    Integer :: m
    REAL(prcn) ::  dfll(nspmx), da(2)

    REAL(prcn) ::  xl(ndim), phil(nspmx),philo(nspmx),philo2(nspmx),phili(nspmx)
    REAL(prcn) ::  nlphil(nspmx),onlphil(nspmx)
    INTEGER :: sp, i, j, k, n, l 
    REAL(prcn) :: rad, tempor(nspmx), gradphibnd(3,nspmx), normal(3), sourcell(nspmx), normal_dotaxial, normal_along_axial(1:3), normal_along_perp(1:3)

    INTEGER :: ib, ie, jb, je, kb, ke, onew, ii, jj, kk 

    Integer :: pcell(3)
    integer ::  is(ndim),iii(ndim),io(ndim), io2(ndim)
    REAL(prcn) ::  tempt, areabdy(nbody), gradphi_tan(3), &
         & phisurface_tmp(nbody,nspmx), area_spec1, flux_tmp,&
         & nusselt_avg(4), source_mean(nspmx), tmp_nu_error_array(nerr_steps), flux_body_parallel(nbody,nspmx), flux_body_perp(nbody,nspmx), gradphi_dotaxial, gradphi_axial(1:3)    

    sum_nlterm = zero 
    ALLOCATE(flux_body_components(ndim, nbody, nspmx),&
         & flux_global_parallel(nspmx), flux_global_perp(nspmx))
    flux_body_components = zero 
    flux_global_parallel = zero
    flux_global_perp = zero 
    flux_body_parallel = zero
    flux_body_perp = zero 
    !CALL calc_gradphi

    DO m = 1, nbody !Loop over bodies

       da(1)=4.*pi*(radbdy(m)*dx)**2./float(nbnd)
       areabdy(m) = 4.*pi*(radbdy(m)*dx)**2
       flux_body(m,:) = zero
       flux_body2(m,:) = zero

       DO  10 l=1,nbnd

          rad = zero
          phil(:)=zero
          dfll(:) = zero
          DO n=1,ndim

             xl(n)=xc(m,n)+xs(n,l)*radbdy(m)

             if(xl(n).lt.zero) then 
                is(n) = int(xl(n)-1)
             else 
                is(n) = int(xl(n))
             end if

             rad=rad+(xs(n,l)*radbdy(m))**2.0
             !Print*,'radius  = ', rad, dsqrt(rad)
          ENDDO

          rad = SQRT(rad)
          normal(1:3) = xs(1:3,l)
          !print*,'normap = ', normal, sqrt(normal(1)**2.d0+normal(2)&
          !     &**2.d0+normal(3)**2.d0)
          !pcell(1:3) = is(1:3)

          
          call interpolate_phidata(is,xl, ib&
               &,ie,jb,je,kb,ke,phil,nlphil,onlphil,dfll, 1,onew) 

          sum_nlterm = sum_nlterm + nlphil(1)*da(1)*dx
          do k = 1, onew
             do j = 1, onew
                do i = 1, onew
                   DO sp = 1, nspmx
                      ii = ib+i-1
                      jj = jb+j-1
                      kk = kb+k-1
                      if(ii.lt.1) ii = mxf+ii-1
                      if(ii.gt.mxf-1) ii = ii-mxf +1
                      if(jj.lt.1) jj = my+jj
                      if(jj.gt.my) jj = jj-my
                      if(kk.lt.1) kk = mz+kk
                      if(kk.gt.mz) kk = kk-mz 

                      gradphisten(i,j,k,:,sp) = gradphi(ii,jj,kk,:)

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


             !flux_body_components(1:ndim, m, sp) = flux_body_components(1:ndim, m, sp) + gradphibnd(1:ndim, sp)*normal(1:ndim)

             normal_dotaxial =  array_dot_product(normal(1:3), axial_direction(1:3))

             normal_along_axial(1:3) = normal_dotaxial*axial_direction(1:3)

             normal_along_perp(1:3) = normal(1:3) - normal_along_axial(1:3)


             gradphi_dotaxial = array_dot_product(gradphibnd(1:3,sp), axial_direction(1:3))

             gradphi_axial(1:3) = gradphi_dotaxial*axial_direction(1:3)
             gradphi_tan(1:3) = gradphibnd(1:3,sp) - gradphi_dotaxial*axial_direction(1:3)

             flux_body_parallel(m,sp) = flux_body_parallel(m,sp) + array_dot_product(gradphi_axial(1:3), normal_along_axial(1:3))


             flux_body_perp(m,sp) = flux_body_perp(m,sp) + array_dot_product(gradphi_tan(1:3), normal_along_perp(1:3))

          end do

          !RINT*,'flux = ', flux_body(m,1)


10     ENDDO    !close loop over all boundary points

       flux_global_parallel(:) = flux_global_parallel(:) + flux_body_parallel(m,:)*da(1)
       flux_global_perp(:) = flux_global_perp(:) + flux_body_perp(m,:)*da(1)
       !flux_global_
    END DO !end the body loop

    nusselt_para(imis) = -(flux_global_parallel(1)*dchar*dchar)&
         &/(6.d0*voldom*maxvolfrac*(phisurf-phistream))
    nusselt_perp(imis) = -(flux_global_perp(1)*dchar*dchar)&
         &/(6.d0*voldom*maxvolfrac*(phisurf-phistream))
    nusselt_all_mis(imis) = -((flux_global_perp(1)&
         &+flux_global_parallel(1))*dchar*dchar)/(6.d0*voldom&
         &*maxvolfrac*(phisurf-phistream)) 


    !WRITE(*,'(A,3(2x,g17.8))') 'FLUX PARALLEL AND PERP = ', flux_global_parallel, flux_global_perp

    !open(unit=200, file = TRIM(RUN_NAME)//'_POST_nusselt.dat', form="formatted", status = "unknown")
    !write(200, '(10(2x,g17.8))') nusselt_parallel, nusselt_perp, nusselt_total

  END subroutine calc_flux

  REAL(prcn) FUNCTION  MAGVEC(vecin)
    IMPLICIT NONE
    INTEGER :: idim, dimvec
    REAL(prcn) , DIMENSION(:), INTENT(IN) :: vecin

    magvec = zero 
    dimvec = size(vecin,1)

    do idim = 1, dimvec
       magvec = magvec + vecin(idim)**two
    end do
    magvec = dsqrt(magvec)
    RETURN
  end FUNCTION MAGVEC


  SUBROUTINE calc_diffn 
    USE nlarrays, ONLY : ff1=>uf1, ff2=>uf2, ff3=>uf3,fr1&
         &=>ur1,fr2=>ur2,fr3=>ur3, ff4 => uf11, ff5 => uf22, ff6 => uf33
    IMPLICIT NONE 
    INTEGER :: sp, i, j, k, ioffset
    Real(prcn) :: sumdiff1, sumdiff2, sumdiff3
    sumdiff1 = zero
    sumdiff2 = zero
    sumdiff3 = zero


    DO sp =1, nspmx ! Initialize the loop for all the species (nspmx)
       DO i=1,mxf !Was 2,mxf previously for reasons unknown RG 02/28/06
          ioffset=foffset+i
          DO k=1,mz
             DO j=1,my2
                if(xperiodic) then
                   if(ioffset.eq.1) THEN
                      ff1(j,k)=(1./dx2)*(phif(mxf-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+ phif(ioffset+1,j,k,sp)) 
                   else if(ioffset.eq.mxf) THEN
                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(2,j,k,sp))

                   else

                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(ioffset+1,j,k,sp)) 
                   endif
	        else
                   ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.&
                        &*phif(ioffset,j,k,sp)+phif(ioffset+1,j,k,sp)) 
                endif


                ff1(j,k)=ff1(j,k)*gamma(sp)

                ff5(j,k)= gamma(sp)*wy(j)*wy(j)*phif(ioffset,j,k,sp)

                ff6(j,k) = gamma(sp)*wz(k)*wz(k)*phif(ioffset,j,k,sp)

                !end of diffusion terms calculation 
                !beginnig of grad phi tems calculation
		if(i.le.mx1)then
                   if(i.eq.1) then
                      ff2(j,k) = (phif(2,j,k,sp)-phif(mxf-1,j,k,sp))/(two*dx) 

                   else if(i.eq.(mxf))then
                      ff2(j,k) = (phif(2,j,k,sp)-phif(i-1,j,k,sp))/(two*dx) 
                   ELSE

                      ff2(j,k) = (phif(i+1,j,k,sp)-phif(i-1,j,k,sp))/(two*dx)
                   end if


                   ff3(j,k)=phif(i,j,k,sp)*wy(j)  !starts at foffset+1
                   ff4(j,k)=phif(i,j,k,sp)*wz(k)
                ENDIF
             ENDDO
          ENDDO

          CALL ff2cr(ff1(:,:),diffn(i,:,:,1))
          CALL ff2cr(ff5(:,:),diffn(i,:,:,2))
          CALL ff2cr(ff6(:,:),diffn(i,:,:,3))

          if(i.le.mx1) then  
             CALL ff2cr(ff2(:,:),gradphi(i,:,:,1))
             CALL ff2cr(ff3(:,:),gradphi(i,:,:,2))
             CALL ff2cr(ff4(:,:),gradphi(i,:,:,3))
          endif
       END DO
    ENDDO
    !write(90,'(12(2x,e20.10))') t, sumdiff1/(mxf*my*mz), sumdiff2/(mxf&
    !     &*my*mz), sumdiff3/(mxf*my*mz)
  END SUBROUTINE calc_diffn


  SUBROUTINE calc_gradphi
    USE nlarrays, ONLY : ff1=>uf1, ff2=>uf2, ff3=>uf3, ff4=>uf11,fr1&
         &=>ur1,fr2=>ur2,fr3=>ur3 
    IMPLICIT NONE 
    INTEGER :: sp, i, j, k, ioffset


    DO sp =1, nspmx ! Initialize the loop for all the species (nspmx)
       DO i=1,mxf !Was 2,mxf previously for reasons unknown RG 02/28/06
          ioffset=foffset+i
          DO k=1,mz
             DO j=1,my2
                !beginnig of grad phi tems calculation
		if(i.le.mx1)then
                   if(i.eq.1) then
                      ff2(j,k) = (phif(2,j,k,sp)-phif(mxf-1,j,k,sp))/(two*dx) 

                   else if(i.eq.(mxf))then
                      ff2(j,k) = (phif(2,j,k,sp)-phif(i-1,j,k,sp))/(two*dx) 
                   ELSE

                      ff2(j,k) = (phif(i+1,j,k,sp)-phif(i-1,j,k,sp))/(two*dx)
                   end if


                   ff3(j,k)=phif(i,j,k,sp)*wy(j)  !starts at foffset+1
                   ff4(j,k)=phif(i,j,k,sp)*wz(k)
                ENDIF
             ENDDO
          ENDDO

          if(i.le.mx1) then  
             CALL ff2cr(ff2(:,:),gradphi(i,:,:,1))
             CALL ff2cr(ff3(:,:),gradphi(i,:,:,2))
             CALL ff2cr(ff4(:,:),gradphi(i,:,:,3))
          endif
       END DO
    ENDDO
  END SUBROUTINE calc_gradphi
end MODULE scalar_post
