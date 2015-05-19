SUBROUTINE histogram(u,wt,n,nddu,bldd,brdd,hist)
  !nddu: number of bins for pdf formation
  !bldd: lefh hand side limit .... output
  !brdd: right side limie.... otuput
  USE randomno
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n , nddu
  REAL(prcn), INTENT(in), DIMENSION(:) :: u, wt
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: u_t

  REAL(prcn),  INTENT(out) :: bldd, brdd

  REAL(prcn),  INTENT(out), DIMENSION(:) ::  hist(100)
  REAL(prcn) ::  xdiff

  REAL(prcn) :: vlmt, ave, adev, sdev, var, skew, curt

  INTEGER :: i, ibin
  bldd= 1.e25
  brdd=-1.e25
  ALLOCATE(u_t(n))
  CALL moment1(4, n, n, wt, u, ave,adev,sdev,var,skew&
       &,curt)
  PRINT*,'ave,var,sdev, skew, curt=', ave, var,sdev, skew,curt
  WRITE(*,*)'number of bins in hist..',nddu

  DO i=1,n
     u_t(i)=(u(i)-ave)/sdev
  ENDDO

  CALL moment1(4, n, n, wt, u_t, ave,adev,sdev,var,skew&
       &,curt)
  PRINT*,'ave,var, sdev,skew, curt=', ave, var,sdev, skew,curt
  DO i=1,nddu
     hist(i)=0.0  
  ENDDO
  bldd = MIN(MINVAL(u_t(:)), bldd)
  brdd = MAX(MAXVAL(u_t(:)), brdd)
!!$    DO i=1,n
!!$       bldd = amin1(u(i),bldd)
!!$       brdd = amax1(u(i),brdd)
!!$    ENDDO

  xdiff  = (brdd - bldd)/float(nddu-1)

  DO i=1,n
     ibin = (u_t(i) - bldd)/xdiff + 1
     hist(ibin)=hist(ibin) + wt(i)/xdiff
  ENDDO

  DEALLOCATE(u_t)

END SUBROUTINE histogram

 SUBROUTINE plothist(hist,lb,ub,nbins,iunit,t,tref)
   IMPLICIT NONE
   REAL(prcn), INTENT(in), DIMENSION(:) ::hist

   INTEGER, INTENT(in) ::iunit,nbins
   REAL(prcn) ::  sum_x,lb,ub,dx, t, tref, tmp, tmp2
   INTEGER(prcn) :: i

   sum_x = lb
   dx= (ub-lb)/(float(nbins)-1)
   tmp = one/sqrt(twopi)
   WRITE(iunit,*)'Zone t="',t/tref,'"'
   DO i=1,nbins
      tmp2  = sum_x+dx/2
      WRITE(iunit,*)tmp2,hist(i)!, tmp*exp(-(tmp2*tmp2)/two)
      sum_x = sum_x + dx
   ENDDO

   RETURN
 END SUBROUTINE plothist

