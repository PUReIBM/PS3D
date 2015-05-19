MODULE test_flux
  USE PRECISION 
  USE constants 
  USE fftw_interface
  USE scalar_data 
  USE nr
  IMPLICIT NONE 
  
  REAL(prcn), PRIVATE :: phireal(mx,my,mz), flux_anal, Nu3_an(nbmx&
       &,nsmx,nspmx), Nu3_anint(nbmx,nsmx,nspmx),phi_anint(nbmx,nsmx&
       &,nspmx,2),phi_anal(nbmx,nsmx&
       &,nspmx,2), alpha

  REAL(prcn) :: Nu3_anint2(nbmx,nsmx,nspmx)
CONTAINS 
  SUBROUTINE calc_anal
    IMPLICIT NONE 
    REAL(prcn) :: dist, time, tmp
    INTEGER :: i,j,k, m
    m=1 
    time  = t! dt*nstep
    !time = 20.d0
    PRINT*,'dt=', dt, time,t 
    alpha = gamma(1)
    DO i = 1,mx
       DO j =1, my
          DO k = 1, mz
             dist=(xc(m,1)+foffset-float(i))**2 +(xc(m,2)-float(j))&
                  &**2+(xc(m,3)-float(k))**2
             dist=SQRT(dist)
             IF((dist.GE.bndrad))THEN
                CALL get_tempanal(dist, time, phireal(i,j,k))
                !tmp = ((dist-bndrad)*dx)/(two*SQRT(alpha*time))
                !phireal(i,j,k) = ((phisurf*bndrad)/(dist))&
                !     &*erfcc(tmp)
                !Print*,'xc =',xc(1,1), foffset
                !print*,'dist = ', dist, bndrad
                !print*,'phi=',phireal(i,j,k), erfcc(tmp), tmp
             ELSE 
                phireal(i,j,k) = phisurf
             END IF
          END DO
       END DO
    END DO
    Print*,'bndrad*dx = ', bndrad*dx
    !read(*,*)
    flux_anal = -(phisurf/(bndrad*dx))*(one+((bndrad*dx)/(SQRT(pi*alpha&
         &*time))))
    
    CALL write_anal
    CALL calc_flux
    
  END SUBROUTINE calc_anal
  SUBROUTINE get_tempanal(dist, time, sol)
    IMPLICIT NONE 
    REAL(prcn), INTENT(in) :: dist, time
    REAL(prcn), INTENT(out) :: sol 
    REAL(prcn) :: tmp 
    tmp = ((dist-bndrad)*dx)/(two*SQRT(alpha*time))
    sol = ((phisurf*bndrad)/(dist))&
         &*erfcc(tmp)
  END SUBROUTINE get_tempanal
  
    
  SUBROUTINE calc_flux
    USE errormesgs
    USE general_funcs
    USE dependent_functions
    USE Interpolation
    USE initialize_flo
    IMPLICIT NONE 
    
    REAL(prcn) :: phir(mxf,my,mz,nspmx), phiranal(mxf,my,mz,nspmx)
    INTEGER :: m , l , n , isp , i, j ,k , sp, iii
    INTEGER ::  is(ndim),ii(ndim),io(ndim),io2(ndim)
    REAL(prcn) :: rad, a1(2,2,2),a12(2,2,2), tnorm,  snorm(ndim), rad2
    REAL(prcn) :: philo(nspmx), phili(nspmx), philo2(nspmx), philo_an(nspmx)&
         &, philo2_an(nspmx), philo_anint(nspmx), philo2_anint(nspmx)
    REAL(prcn) ::  t1(ndim),t1o(ndim),t1i(ndim),t1o2(ndim)
    REAL(prcn) ::  xl(ndim),xlo(ndim),xli(ndim),xlo2(ndim)
    
    INTEGER :: ib, ie, jb, je, kb, ke, pcell(3), onew, unitno
    REAL(prcn) ::  dphiint(nspmx,ndim)
    
    REAL(prcn) :: tmp1, tmp2
    DO i=foffset+1,mxf+foffset
       iii = i-foffset 
       CALL ff2cr(phif(i,:,:,1),phir(iii,:,:,1))
       phiranal(iii,:,:,1) = phireal(i,:,:)
    END DO
        
    DO m=1,nbody            ! loop over bodies
       DO l=1,nbnd
          rad = zero
          
          DO n=1,ndim
             dphiint(:,n)=zero
          ENDDO
          DO n=1,ndim
             
             xl(n)=xc(m,n)+xs(n,l)
             
             is(n)=INT(xl(n))
             t1(n)=xl(n)-is(n)
             IF (t1(n).GT.1.) WRITE(*,*) 't1 overflow'
             IF (t1(n).LT.0.) WRITE(*,*) 't1 underflow'
             
             !if(lc.eq.nstep)
             rad=rad+(xs(n,l)*xs(n,l))
             !Print*,'radius  = ', rad, dsqrt(rad)
          ENDDO

          !if(lc.eq.nstep)then
          rad =dsqrt(rad)
          
          DO n=1,ndim
             snorm(n)=xs(n,l)/rad
          END DO
          
          DO k=0,1
             DO j=0,1
                DO i=0,1
                   
                   !-----------------------------------------------------
                   !     Compute the interpolation functions for u and p 
                   !-----------------------------------------------------
                   
                   a1(i+1,j+1,k+1)=((1-i)+(2*i-1)*t1(1))* ((1-j)+(2&
                        &*j-1)*t1(2))*((1-k)+(2*k-1)*t1(3)) 
                ENDDO
             ENDDO
          ENDDO
          
          DO k=0,1
             DO j=0,1
                DO i=0,1
                   DO sp=1,nspmx
                      DO n=1,ndim
                         dphiint(sp,n)= dphiint(sp,n)+a1(i+1,j+1,k&
                              &+1)* dphir(is(1)+i,is(2)+j,is(3)+k,n&
                              &,sp) 
                      ENDDO
                      
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          tnorm = zero
          DO sp=1,nspmx       !check for > 1 scalars
             DO n=1,ndim
                tnorm = tnorm + dphiint(sp,n)*snorm(n)
             ENDDO
             Nu2(m,l,sp)=(tnorm)
             
             !PRINT*, dphiint(1,1),dphiint(1,2), tnorm
          ENDDO
       END DO
       DO l=1,nrpr
        
          rad = zero
          rad2 = zero
          DO sp=1,nspmx
             philo(sp)=zero
             philo2(sp)=zero
             philo_anint(sp)=zero
             philo2_anint(sp)=zero
             
          ENDDO
          
          DO n=1,ndim
             
             !     location of internal points

             xli(n)=xc(m,n)+xi(n,l)

             ii(n)=INT(xli(n))
             
             !     location of external points
             
             xlo(n)=xc(m,n)+xo(n,l)
             xlo2(n) = xc(m,n)+xo2(n,l)
             
             io(n)=INT(xlo(n))
             io2(n)=INT(xlo2(n))
             
             t1i(n)=xli(n)-ii(n)

             t1o(n)=xlo(n)-io(n)
             t1o2(n)=xlo2(n)-io2(n)
             rad=rad+(xo(n,l)*xo(n,l))
             rad2=rad2+((xo2(n,l)**2.))

          ENDDO
          


          
          !Print*,'xlo = ', xlo(1), xlo(2),xlo(3)
          !Print*, 'xc = ', xc(m,1), xc(m,2), xc(m,3)
          !Print*, 'tlo = ', t1o(1), t1o(2)
          !Read(*,*)
          rad=dsqrt(rad)
          rad2=dsqrt(rad2)
          !print*,'rado', rad,' rado2 = ', (rad2)
          DO k=0,1
             DO j=0,1
                DO i=0,1

                   a1(i+1,j+1,k+1)=((1-i)+(2*i-1)*t1o(1))*((1-j)+(2&
                        &*j-1)*t1o(2))*  ((1-k)+(2*k-1)*t1o(3)) 

                   a12(i+1,j+1,k+1)=((1-i)+(2*i-1)*t1o2(1))*((1-j)+(2&
                        &*j-1)*t1o2(2))*((1-k)+(2*k-1)*t1o2(3))
              
                   DO sp=1,nspmx
                      philo(sp)=philo(sp)+a1(i+1,j+1,k+1)*&
                           & phir(io(1)+i,io(2)+j,io(3)+k,sp) 
                      philo2(sp)=philo2(sp)+a12(i+1,j+1,k+1)*&
                           & phir(io2(1)+i,io2(2)+j,io2(3)+k,sp)
                      
                      philo_anint(sp)=philo_anint(sp)+a1(i+1,j+1,k+1)*&
                           & phiranal(io(1)+i,io(2)+j,io(3)+k,sp) 
                      philo2_anint(sp)=philo2_anint(sp)+a12(i+1,j+1,k+1)*&
                           & phiranal(io2(1)+i,io2(2)+j,io2(3)+k,sp) 
                   ENDDO
                   !Print*,i,j,k,a1(i+1,j+1,k+1), phiranal(io(1)+i&
                   !     &,io(2)+j,io(3)+k,1), philo_anint(1)      
                ENDDO
             ENDDO
          ENDDO
          pcell = io
          tmp1  = philo(1)
          tmp2  = philo2(1)
          !tmp1 = philo_anint(1)
          !tmp2 = philo2_anint(1)
          
          CALL set_interpolation_scheme(2) !for LNI
          !Print*,'pcell =',pcell(1), pcell(2),'philo M1='&
          !     &,philo_anint(1)
          !Print*,'interpolation scheme=', interp_scheme,' order = ',order
          CALL set_interpolation_stencil(pcell,ib,ie,jb,je,kb,ke&
               &,interp_scheme, onew) 
          !Print*,'onew = ', onew, size(gstencil,1), size(gstencil,2)
          DO  k = 1, onew 
             DO j = 1,onew
                DO i = 1,onew
                   gstencil(i,j,k,1) = ib+(i-1)
                   gstencil(i,j,k,2) = jb+(j-1)
                   gstencil(i,j,k,3) = kb+(k-1)
                ENDDO
             ENDDO
          ENDDO
          
          vstencil(1:onew,1:onew,1:onew,1) = phir(ib:ie,jb:je,kb:ke,1)
          
          CALL interpolator(gstencil,vstencil(1:onew,1:onew,1:onew,1)&
               &,xlo,philo(1),onew, interp_scheme,weightp)
      
          !
          !Now do it for the second layer of reversal points.
          !
          CALL set_interpolation_scheme(2) !
          
          pcell = io2
          CALL set_interpolation_stencil(pcell,ib,ie,jb,je,kb,ke&
               &,interp_scheme, onew) 
          !Print*,'onew = ', onew, size(gstencil,1), size(gstencil,2)
          DO  k = 1, onew 
             DO j = 1,onew
                DO i = 1,onew
                   gstencil(i,j,k,1) = ib+(i-1)
                   gstencil(i,j,k,2) = jb+(j-1)
                   gstencil(i,j,k,3) = kb+(k-1)
                ENDDO
             ENDDO
          ENDDO
          
          vstencil(1:onew,1:onew,1:onew,1) = phir(ib:ie,jb:je,kb:ke,1)
          
          CALL interpolator(gstencil,vstencil(1:onew,1:onew,1:onew,1)&
               &,xlo2,philo2(1),onew, interp_scheme,weightp)
          
          DO sp=1,nspmx
             CALL get_tempanal(rad, t, philo_an(sp))
             CALL get_tempanal(rad2, t, philo2_an(sp))
             
             phi_anint(m,l,sp,1) = philo_anint(sp)
             phi_anint(m,l,sp,2) = philo2_anint(sp)
             
             !phi_anal(m,l,sp,1) = philo(sp)
             !phi_anal(m,l,sp,2) = philo2(sp)
             phi_anal(m,l,sp,1) = philo_an(sp)
             phi_anal(m,l,sp,2) = philo2_an(sp)
             
             !Print*,'phi=', phi_anint(m,l,sp,1), phi_anal(m,l,sp,1)
             Nu1(m,l,sp)=(-phisurf+philo(sp))/((-bndrad&
                  &+rad)*dx)
             Nu3(m,l,sp)=(-three*phisurf+four*philo(sp)-philo2(sp))&
                  &/(two*(-bndrad+rad)*dx)
             Nu3_an(m,l,sp)=(-three*phisurf+four*philo_an(sp)-philo2_an(sp))&
                  &/(two*(-bndrad+rad)*dx)
             Nu3_anint(m,l,sp)=(-three*phisurf+four*philo_anint(sp)&
                  &-philo2_anint(sp))/(two*(-bndrad+rad)*dx) 
             
             !Print*,  phisurf,phi_anint(m,l,sp,1),phi_anint(m,l,sp,2)&
             !     &, Nu3_an(m,l,sp), (two*(-bndrad+rad)*dx)
             Nu3_anint2(m,l,sp)=(-three*phisurf+four&
               &*tmp1-tmp2)/(two*(-bndrad+rad)&
               &*dx)
          ENDDO
          !Nu2 is mix of FD and Fourier differentiation called FM3
          !Nu1 is central-diff FD called FM1 
          !Nu3 is one sided second order FD scheme called FM2
!!$          Print*,'FIrst rversal pt'
!!$          Print*,'philo M1=',tmp1, philo_anint(1),philo_an(1)
!!$          Print*,'philo M2=',tmp2 , philo2(1), philo_an(1), &
!!$               &(philo_anint(1)-philo_an(1))/philo_an(1)
!!$          Print*,'second rversal pt'
!!$          Print*,'philo2 M1=',tmp2, philo2_an(1), (tmp2&
!!$               &-philo2_an(1))/philo2_an(1),(-three*phisurf+four&
!!$               &*philo_anint(1)-tmp2)/(two*(-bndrad+rad)&
!!$               &*dx), Nu3_an(m,l,1)   
!!$          Print*,'philo2 M2=',philo2_anint(1), philo2_an(1), &
!!$               &(philo2_anint(1)-philo2_an(1))/philo2_an(1),(-three&
!!$               &*phisurf+four*philo_anint(1)-philo2_anint(1))/(two&
!!$               &*(-bndrad+rad)*dx)   
       END DO
    END DO!nbody
    
!!$    unitno = getnewunit(minunitno,maxunitno)
!!$    
!!$    IF (unitno.LT.0) CALL printerror("newunit","ounit")
!!$
!!$    OPEN(unit=unitno,file='lpivslni.dat',form='formatted',status='repl&
!!$         &ace')
!!$
!!$    DO m=1,nbody
!!$       WRITE(unitno,*)'Zone'
!!$       DO isp = 1,nspmx
!!$          DO l=1,nrpr
!!$             IF (xo(2,l).GE.zero.AND.xo(3,l).EQ.zero) THEN 
!!$                IF(xo(1,l).GE.zero) THEN 
!!$                   WRITE(unitno,21)180.-((180.*ATAN(xo(2,l)/xo(1,l)))&
!!$                        &/pi),Nu3_anint(m,l,isp),Nu3_anint2(m,l,isp),&
!!$                        & flux_anal, tmp2, philo2_anint(1),&
!!$                        & philo2_an(1),philo_anint(1),philo_an(1) 
!!$                ELSE
!!$                   WRITE(unitno,21)((180. *ATAN(-xo(2,l)/xo(1,l)))&
!!$                        &/pi),Nu3_anint(m,l,isp), Nu3_anint2(m,l,isp)&
!!$                        &,flux_anal, tmp2, philo2_anint(1),&
!!$                        & philo2_an(1),philo_anint(1),philo_an(1)
!!$                ENDIF
!!$             ENDIF
!!$          ENDDO
!!$       ENDDO
!!$    END DO
21  FORMAT(10(1xe17.4))
!    CLOSE(unitno,status= 'keep')
    PRINT*,'write out the Nu nos'
    CALL write_flux 
    call calc_error
    !call write_phianint
  END SUBROUTINE calc_flux
  
  SUBROUTINE calc_error
    USE errormesgs
    USE general_funcs
    IMPLICIT NONE 
    REAL(prcn) :: err_nu1, err_nu2, err_nu3, err_an, err_anint, err_anint2
    Integer :: m,l,sp, unitno
    err_nu1 = zero
    err_nu2 = zero
    err_nu3 = zero
    err_an = zero
    err_anint = zero
    err_anint2 = zero
    
    !Nu2 is mix of FD and Fourier differentiation called FM3
    !Nu1 is central-diff FD called FM1 
    !Nu3 is one sided second order FD scheme called FM2
    do m = 1,nbody
       do l = 1,nrpr
          do sp = 1, nspmx
             !Print*,'nus = ',Nu1(m,l,sp), flux_anal
             err_nu1 = err_nu1+abs((Nu1(m,l,sp)-flux_anal)/flux_anal)!**2.0
             err_nu2 = err_nu2+abs((Nu2(m,l,sp)-flux_anal)/flux_anal)!**2.0
             err_nu3 = err_nu3+abs((Nu3(m,l,sp)-flux_anal)/flux_anal)!**2.0
             err_an = err_an + abs((Nu3_an(m,l,sp)-flux_anal)/flux_anal)!**2.
             err_anint = err_anint + abs((Nu3_anint(m,l,sp)-flux_anal)&
                  &/flux_anal)!**2. 
             !Print*, 'err = ',abs((Nu3_an(m,l,sp)-flux_anal)&
             !     &/flux_anal), l, flux_anal
             err_anint2 = err_anint2 + abs((Nu3_anint2(m,l,sp)&
                  &-flux_anal)/flux_anal)!**2.
             !Print*, err_nu1, err_nu2, err_nu3
          end do
       end do
    end do
    err_nu1 = err_nu1/nrpr
    err_nu2 = err_nu2/nrpr
    err_nu3 = err_nu3/nrpr
    err_an = err_an/nrpr
    err_anint = err_anint/nrpr
    err_anint2 = err_anint2/nrpr
    unitno = getnewunit(minunitno,maxunitno)
    
    IF (unitno.LT.0) CALL printerror("newunit","ounit")
    
    OPEN(unit=unitno,file='err_flux.dat',form='formatted',status='replace')
    write(unitno, 21) dx, err_nu1, err_nu3, err_nu2, err_anint2,&
         & err_an, err_anint 
    close(unitno, status='keep')
    21 format(10(1x,e17.4))
  end SUBROUTINE calc_error

  SUBROUTINE write_flux
    USE errormesgs
    USE general_funcs
    IMPLICIT NONE 
    INTEGER :: i,j,k,m,l,isp, unitno
    unitno = getnewunit(minunitno,maxunitno)
    
    IF (unitno.LT.0) CALL printerror("newunit","ounit")

    OPEN(unit=unitno,file='thetavsNU1.dat',form='formatted',status='re&
         &place')
    
    !Nu2 is mix of FD and Fourier differentiation called FM3
    !Nu1 is central-diff FD called FM1 
    !Nu3 is one sided second order FD scheme called FM2
    DO m=1,nbody
       WRITE(unitno,*)'Zone'
       DO isp = 1,nspmx
          DO l=1,nrpr
             IF (xo(2,l).GE.zero.AND.xo(3,l).EQ.zero) THEN 
                IF(xo(1,l).GE.zero) THEN 
                   WRITE(unitno,21)180.-((180.*ATAN(xo(2,l)/xo(1,l)))/pi)&
                        &,-Nu1(m,l,isp)*two*bndrad*dx,-Nu3(m,l,isp)*two*bndrad*dx,-Nu2(m,l,isp)*two*bndrad*dx,- Nu3_anint2(m,l,isp)*two*bndrad*dx&
                        &,-Nu3_an(m,l,isp)*two*bndrad*dx,-Nu3_anint(m,l,isp)*two*bndrad*dx, -flux_anal*two*bndrad*dx
                ELSE
                   WRITE(unitno,21)((180. *ATAN(-xo(2,l)/xo(1,l)))/pi)&
                        &,-Nu1(m,l,isp)*two*bndrad*dx ,-Nu3(m,l,isp)*two*bndrad*dx, -Nu2(m,l,isp)*two*bndrad*dx, -Nu3_anint2(m,l,isp)*two*bndrad*dx&
                        &,-Nu3_an(m,l,isp)*two*bndrad*dx, -Nu3_anint(m,l,isp)*two*bndrad*dx,-flux_anal*two*bndrad*dx
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    END DO
21  FORMAT(10(1xe17.4))
    CLOSE(unitno,status= 'keep')
  END SUBROUTINE write_flux

  
  SUBROUTINE  write_phianint 
    USE errormesgs
    USE general_funcs
    IMPLICIT NONE 
    INTEGER :: i,j,k,m,l,isp, unitno
    unitno = getnewunit(minunitno,maxunitno)
    
    IF (unitno.LT.0) CALL printerror("newunit","ounit")
    
    OPEN(unit=unitno,file='phianint.dat',form='formatted',status='replace')
    DO m=1,nbody
       WRITE(unitno,*)'Zone'
       DO isp = 1,nspmx
          DO l=1,nrpr
             IF (xo(2,l).GE.zero.AND.xo(3,l).EQ.zero) THEN 
                IF(xo(1,l).GE.zero) THEN 
                   WRITE(unitno,21)180.-((180.*ATAN(xo(2,l)/xo(1,l)))/pi)&
                        &,phi_anint(m,l,isp,1), phi_anint(m,l,isp,2)&
                        &,phi_anal(m,l,isp,1), phi_anal(m,l,isp,2) 
                ELSE
                   WRITE(unitno,21)((180.*ATAN(-xo(2,l)/xo(1,l)))/pi)&
                        &,phi_anint(m,l,isp,1), phi_anint(m,l,isp,2)&
                        &,phi_anal(m,l,isp,1), phi_anal(m,l,isp,2) 
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    END DO
21  FORMAT(10(1xe17.4))
    CLOSE(unitno,status= 'keep')
    
  END SUBROUTINE write_phianint


  
  SUBROUTINE write_anal
    USE errormesgs
    USE general_funcs
    USE nlarrays , ONLY : ur1,ur2,ur3, phir1=>ur11
    USE fftw_interface
    IMPLICIT NONE 
    REAL(prcn) :: dist
    INTEGER :: i,j,k, unitno, m,ii
    unitno = getnewunit(minunitno,maxunitno)
    
    IF (unitno.LT.0) CALL printerror("newunit","ounit")

    OPEN(unit=unitno,file='phireal.dat',form="formatted"  &
         ,status="replace")
    WRITE(unitno,*)'VARIABLES= ',' "X" ',' "Y" ',' "phi_anal" ',' "phi&
         &_num" ',' "err" '
    WRITE(unitno,*)'ZONE F=POINT, I=', mx,  ', J=', my
    k = mz/2
    DO j = 1, my 
       DO i = 1,mx 
          m=1
          CALL ff2cr(phif(i,:,:,1), phir1)
          dist=(xc(m,1)+foffset-float(i))**2 +(xc(m,2)-float(j))&
               &**2+(xc(m,3)-float(k))**2
          dist=SQRT(dist)
          IF((dist.GE.bndrad.AND.phireal(i,j,k).GT.zero))THEN
             WRITE(unitno, 21)(i-1)*dx,(j-1)*dy, phireal(i,j,k),&
                  & phir1(j,k), (ABS(phireal(i,j,k)-phir1(j,k)))!/phireal(i,j,k)
          ELSE
             WRITE(unitno,21)(i-1)*dx,(j-1)*dy, zero, zero, zero
          END IF
          
       END DO
    END DO
    21 format(10(2xe20.10))
    CLOSE(unitno, status='keep')

    OPEN(unit=unitno,file='dphireal.dat',form="formatted"  &
         ,status="replace")
    WRITE(unitno,*)'VARIABLES= ',' "X" ',' "Y" ',' "dphidx" ',' "dphid&
         &y" ',' "dphidz" ' 
    WRITE(unitno,*)'ZONE F=POINT, I=', mxf,  ', J=', my
    k = mz/2
    DO j = 1, my 
       DO i=foffset+1,mxf+foffset
          ii = i-foffset 
          m=1
          WRITE(unitno,*)(i-1)*dx,(j-1)*dy, dphir(ii,j,k,1,1)&
               &,dphir(ii,j,k,2,1), dphir(ii,j,k,3,1) 
       END DO
    END DO
    CLOSE(unitno, status='keep')
  END SUBROUTINE write_anal

  
END MODULE test_flux


