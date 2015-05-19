MODULE nl_allflow_noncons
  USE precision 
  USE constants 
  USE global_data
  USE fftw_interface
  Use nlarrays, dotur11=>ur12, dotur12=>ur13, dotur13=>ur23, dufdx=>uf11, dufdy=>uf22, dufdz=>uf33, ff1=>uf12, ff2=>uf13, ff3=>uf23
  USE bcsetarrays, dur=>fr

  IMPLICIT NONE 
  Private
!!$  COMPLEX(prcn) ::  uf1(my2,mz),uf2(my2,mz),uf3(my2,mz)
!!$  COMPLEX(prcn) ::  uf11(my2,mz),uf22(my2,mz),uf33(my2,mz)
!!$  COMPLEX(prcn) ::  uf12(my2,mz),uf13(my2,mz),uf23(my2,mz)
!!$  COMPLEX(prcn) ::  tmp
!!$  REAL(prcn) :: ur1(my,mz),ur2(my,mz),ur3(my,mz)
!!$  REAL(prcn) :: ur11(my,mz),ur22(my,mz),ur33(my,mz)
!!$  REAL(prcn) :: ur12(my,mz),ur13(my,mz),ur23(my,mz)
  INTEGER :: mx2
  Public :: form_nl_noncons
CONTAINS 
  SUBROUTINE form_nl_noncons!(ubc,pbc)

    Use nlmainarrays, Only : ubcp,pbcp, nlbcp
    IMPLICIT NONE 
    INTEGER :: i,j,k, i1, i2, i3, idim

    REAL(PRCN) ::  u_max, v_max, w_max

    REAL(PRCN) ::  umean_tmp(ndim)

    REAL(PRCN) :: umax_tmp


    External fplus, fminus
!!$    real(prcn), Intent(out) ::  ubc(:,:,:,:)
!!$    real(prcn), intent(out) ::  pbc(:,:,:)


    DTNEW = LARGE_NUMBER
    u_max = SMALL_NUMBER
    v_max = SMALL_NUMBER
    w_max = SMALL_NUMBER
    umax_tmp = LARGE_NUMBER

    mx2 = mx+1
    do idim = 1, ndim 
       do i = 1, mx1
          DO k = 1,mz
             DO j = 1, my2
                if(idim.eq.1) then 
                   uf1(j,k) = u(i,j,k,1) !+ umean(1)

                   uf2(j,k) = u(i,j,k,2)! + umean(2)

                   uf3(j,k) = u(i,j,k,3)! + umean(3)

                   if(i.eq.1) then
                      ff1(j,k) = (p(i+1,j,k)-p(mxf-1,j,k))/(2.*dx) 

                   else if(i.eq.mx1)then
                      ff1(j,k) = (p(1,j,k)-p(i-1,j,k))/(2.*dx) 
                   ELSE

                      ff1(j,k) = (p(i+1,j,k)-p(i-1,j,k))/(two*dx)
                   end if

                   ff2(j,k)=p(i,j,k)*wy(j)  !starts at foffset+1
                   ff3(j,k)=p(i,j,k)*wz(k)
                end if

                !dufdx(j,k) = (u(i+1,j,k,idim) - u(i-1,j,k,idim))*(half/dx)

                dufdy(j,k) = wy(j)*(u(i,j,k,idim))
                dufdz(j,k) = wz(k)*(u(i,j,k,idim))
             END DO
          END DO
          IF(IDIM.EQ.1) THEN 
             CALL ff2cr(uf1,ubcp(i,:,:,1))
             CALL ff2cr(uf2,ubcp(i,:,:,2))
             CALL ff2cr(uf3,ubcp(i,:,:,3))
             CALL ff2cr(ff1(:,:),ppr(i,:,:,1))
             CALL ff2cr(ff2(:,:),ppr(i,:,:,2))
             CALL ff2cr(ff3(:,:),ppr(i,:,:,3))
             !call ff2cr(p(i,:,:),pbcp(i,:,:))
             do k = 1, mz
                do j = 1, my 
                   ubcp(i,j,k,1) = ubcp(i,j,k,1) + umean(1)
                   ubcp(i,j,k,2) = ubcp(i,j,k,2) + umean(2)
                   ubcp(i,j,k,3) = ubcp(i,j,k,3) + umean(3)

                   ppr(i,j,k,1) = ppr(i,j,k,1) + mpg(1)
                   ppr(i,j,k,2) = ppr(i,j,k,2) + mpg(2)
                   ppr(i,j,k,3) = ppr(i,j,k,3) + mpg(3)

                end do
             end do

          end IF
          !CALL ff2cr(dufdx,dur(i,:,:,1))
          CALL ff2cr(dufdy,dur(i,:,:,2))
          CALL ff2cr(dufdz,dur(i,:,:,3))

       end do
       ubcp(mx,:,:,:) = ubcp(1, :,:,:)

       do i = 1, mx1
          DO k = 1,mz
             DO j = 1, my
                IF(IDIM.EQ.1) THEN 
                   if(FLUID_ATIJK(i,j,k))then

                      UMEAN_TMP(1)= UMEAN_TMP(1)+ ubcp(i,j,k,1)!+umean(1)
                      UMEAN_TMP(2)= UMEAN_TMP(2)+ ubcp(i,j,k,2)!+umean(2)
                      UMEAN_TMP(3)= UMEAN_TMP(3)+ ubcp(i,j,k,3)!+umean(3)
                   end if

                   u_max = MAX(u_max, ABS(ubcp(i,j,k,1)))
                   v_max = MAX(v_max, ABS(ubcp(i,j,k,2)))
                   w_max = MAX(w_max, ABS(ubcp(i,j,k,3)))

                end IF

                nlbcp(i,j,k,idim) = (ubcp(i,j,k,2))*dur(i,j,k,2)
                nlbcp(i,j,k,idim) = nlbcp(i,j,k,idim) + (ubcp(i,j,k,3))*dur(i,j,k,3)

                if(ubcp(i,j,k,1).gt.zero) THEN 
                   i1 = i
                   i2 = i-1
                   i3 = i-2
                   if(i2.lt.1) i2 = mxf+i2-1
                   if(i3.lt.1) i3 = mxf+i3-1

                   nlbcp(i,j,k,idim) = nlbcp(i,j,k,idim) + (ubcp(i,j,k,1))*(3.d0*ubcp(i1,j,k,idim)-4.d0*ubcp(i2,j,k,idim)+ubcp(i3,j,k,idim))*(half/dx)
                   !nlbcp(i,j,k,idim) = nlbcp(i,j,k,idim) + (ubcp(i,j,k,1))*(ubcp(i1,j,k,idim)-ubcp(i2,j,k,idim))*(one/dx)
                elseif(ubcp(i,j,k,1).le.zero) THEN 
                   i1 = i
                   i2 = i+1
                   i3 = i+2
                   if(i2.gt.mxf-1) i2 = i2-(mxf-1)
                   if(i3.gt.mxf-1) i3 = i3-(mxf-1)

                   nlbcp(i,j,k,idim) = nlbcp(i,j,k,idim) + (ubcp(i,j,k,1))*(-3.d0*ubcp(i1,j,k,idim)+4.d0*ubcp(i2,j,k,idim)-ubcp(i3,j,k,idim))*(half/dx)
                   !nlbcp(i,j,k,idim) = nlbcp(i,j,k,idim) + (ubcp(i,j,k,1))*(-ubcp(i1,j,k,idim)+ubcp(i2,j,k,idim))*(one/dx)

                end if
                nlbcp(i,j,k,idim) = -nlbcp(i,j,k,idim) !taking nl to the RHS and absorbing the negative sign 
             end DO
          end DO
          call ff2rc(nlbcp(i,:,:,idim), nl(i,:,:,idim))
       end do

       nlbcp(mx,:,:, idim) = nlbcp(1,:,:,idim)
       nl(mx,:,:, idim) = nl(1,:,:,idim)
    end do
    umax_tmp = MAX(u_max,v_max)

    umax_tmp = MAX(umax_tmp,w_max)

    !umax_tmp = MAX(umax_tmp,w_max+umean(3))

    ! if(umax_tmp*dt/dx.gt.cfl) then
    dtnew = cfl*dx/umax_tmp
    !else
    !  dtnew = dt 
    !endif

    PRINT*, "MAX CFL and dtnew", umax_tmp*dt/dx, dt, dtnew
    UMEAN_TMP(:)= UMEAN_TMP(:)/(mx1*my*mz)
    umodmean = DSQRT(umean_tmp(1)**two + umean_tmp(2)**two + umean_tmp(3)**two )

    PRINT*,"UMEAN_TMP = ", UMEAN_TMP, UMODMEAN

!!$    open(1000, file="nl.dat", form="formatted")
!!$    
!!$    write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" ',' "FY" ',' "FZ" '
!!$    
!!$    write(1000,*)'ZONE F=POINT, I=', mxf,  ', J=', my, ', K=', mz

!!$    do k=1,mz
!!$       do j=1,my
!!$          do i=1,mxf
!!$             write(1000,'(10(2x,g12.5))')(real(i)),(real(j)),(real(k)),(nlbcp(i,j,k,1))&
!!$                  &,(nlbcp(i,j,k,2)),(nlbcp(i,j,k,3))
!!$             
!!$          end do
!!$       end do
!!$    end do
!!$    close(1000,status = "keep")


    return 
  end SUBROUTINE form_nl_noncons

  
end MODULE nl_allflow_noncons
