module gauss_seidel
#include "ibm.h"
  ! *********************************************
  USE precision
  USE constants
  USE global_data 
  IMPLICIT NONE


contains

  subroutine gaussseidelp(a,b,c,r,u,n)!,j,k) ! for periodic boundaries and (n=mx-1)
    implicit none
    integer, intent(in) :: n
    real(prcn), Dimension(:), Intent(in)  :: a(n), b(n), c(n), r(n)
    real(prcn), Dimension(:), Intent(inout) :: u(0:n+1)
    real(prcn) , Dimension(:) :: old(0:n+1), firstu(0:n+1)
!!$    INTEGER, OPTIONAL, INTENT(in):: j,k
    integer i
    real(prcn) :: resd, res1, umax, MINTOL, norm2_r, err, tol_gs,norm2_rloc
    LOGICAL ::  ALREADY_PRINT, ITER_EXCEED
!    umax = maxval(ABS(u))
    ! gauss siedel for momentum equations
    
    resd  = one 
    res1 = one 
    iter_gauss = 0
    firstu(:) = u(:) 
    norm2_r = zero
    norm2_rloc = zero
    err = zero 
    ITER_EXCEED = .FALSE.
    ALREADY_PRINT = .FALSE.
    
    do i = 1, n 
       norm2_rloc = norm2_rloc + ( r(i))**two
    end do
    GLOBAL_DOUBLE_SUM(norm2_rloc,norm2_r,1,decomp_group)
    norm2_r = sqrt(norm2_r)
!!$    PRINT*,'norm2_r =', norm2_r
    call  residue(a,b,c,r,u(0:n+1),err,n) 
!!$    if (j.eq.1.and.k.eq.1)then
!!$       PRINT*, 'err = ', err
!!$    end if

!!$ 
!!$    if(resd_r.eq.zero) THEN 
!!$       MINTOL = 1.0E-16
!!$    ELSE
!!$       MINTOL = (1.0E-4)*(umax+SMALL_NUMBER)
!!$    end if
    if(gauss_scal) then 
       tol_gs = 1.0E-14*norm2_r
       
       !WRITE(*,'(A,3(2x,g17.8))') 'norm2_r = ', norm2_r, tol_gs, err
    ELSE
       tol_gs = 1.0E-6*norm2_r
    end if
    
    DO WHILE(err.GT.tol_gs)!1.0e-9 was the last one 
       resd = zero 
       res1 = zero 
10     old(:) = u(:)
       iter_gauss = iter_gauss + 1
 !forward
#if PARALLEL
       do i = 1,n
          u(i) = (r(i) - a(i)*u(i-1) - c(i)*u(i+1))/b(i)
          if(i.eq.1)then
             RSENDRECV(u(i),1,fromproc,1,u(n+1),1,toproc,1,decomp_group,status) 
          else if(i.eq.n)then
             RSENDRECV(u(i),1,toproc,1,u(0),1,fromproc,1,decomp_group,status) 
          end if
       end do
#else
      
       u(1) = (r(1) - a(1)*u(n) - c(1)*u(2))/b(1) ! treat periodic boundaries specially
       do i=2,n-1
          u(i) = (r(i) - a(i)*u(i-1) - c(i)*u(i+1))/b(i)
       enddo
       u(n) = (r(n) - a(n)*u(n-1) - c(n)*u(1))/b(n) ! treat periodic boundaries specially
#endif

       !reverse
#if PARALLEL
       do i=n-1,2
          u(i) = (r(i) - a(i)*u(i-1) - c(i)*u(i+1))/b(i)
       end do
#else
       do i=n-1,2
          u(i) = (r(i) - a(i)*u(i-1) - c(i)*u(i+1))/b(i)
       enddo
       u(1) = (r(1) - a(1)*u(n) - c(1)*u(2))/b(1) ! treat periodic boundaries specially
#endif

       !do i=1, n
       !   resd =  resd + (u(i) - old(i))**2.D0
       !   res1 = res1 +ABS(u(i) - old(i))
       !enddo
       !resd = sqrt(resd)/real(n,prcn)
       !IF(ITER.GE.1) PRINT*,'IN GAUSS', iter, gauss_u, gauss_p, gauss_phi
       ! IF(iter.gt.10000.AND.(.NOT.ALREADY_PRINT)) THEN 

       !PRINT*, 'which gauss =',gauss_u, gauss_p, gauss_phi
!!$          IF(present(j).AND.present(k))then
!!$             PRINT*,  'j = ', j , 'k = ', k
!!$          end IF
!!$       
       !PRINT*, 'residue  = ',err, 1.0E-6*norm2_r, norm2_r, iter
       !PRINT*, ' u = ', u
       !PRINT*, ' firstu = ', firstu
       !PRINT*, ' r = ', r
       !iter = 0
       !ALREADY_PRINT = .TRUE.
       
       !   ITER_EXCEED = .TRUE.
       !end IF
       call  residue(a,b,c,r,u,err,n) 

    end DO
!!$    if(gauss_scal) then 
!!$       
!!$       WRITE(*,'(A,2(2x,e17.8),2x,i10)') 'ERR',tol_gs, err, iter
!!$    end if
!!$    do i = 1,n
!!$       WRITE(*,*)'SOL = ', i, u(i)
!!$    end do
    IF(ITER_EXCEED) PRINT*,'FINAL No. OF ITERS = ', iter_gauss

    !    if(resd.GE.real(1E-6)*(umax+SMALL_NUMBER)) goto 10      

    !if(resd.GE.real(1E-20,prcn)) goto 10      
  End Subroutine gaussseidelp

  Subroutine residue(a,b,c,r,u,err,n)
    Implicit None 
    Integer, intent(in) :: n 
    real(prcn), Dimension(:), Intent(in)  :: a(n), b(n), c(n), r(n)
    real(prcn), Dimension(:), Intent(in) :: u(0:n+1)
    real(prcn), Intent(out) :: err
    Integer :: i 
    real(prcn) :: errloc
    errloc = zero 
    
    do i = 1, n 
#if PARALLEL
       errloc = errloc + (r(i)- (a(i)*u(i-1)+b(i)*u(i)+c(i)*u(i+1)))**two
#else
       If(i.eq.1) then 
          errloc = errloc + (r(i)- (a(i)*u(n)+b(i)*u(i)+c(i)*u(i+1)))**two
       else If(i.eq.n) then 
          errloc = errloc + (r(i)- (a(i)*u(i-1)+b(i)*u(i)+c(i)*u(1)))**two
       else 
          errloc = errloc + (r(i)- (a(i)*u(i-1)+b(i)*u(i)+c(i)*u(i+1)))**two
       end If
#endif
    end do
    GLOBAL_DOUBLE_SUM(errloc,err,1,decomp_group)
!    PRINT*,'errloc =', myid, errloc
    err = sqrt(err)
 !   PRINT*,'err =', myid, err
  !  PARALLEL_FINISH()
   ! STOP
  end Subroutine residue

end module gauss_seidel

    
