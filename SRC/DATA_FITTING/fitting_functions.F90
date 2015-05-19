MODULE fitting_functions
  IMPLICIT NONE 
CONTAINS
  FUNCTION funpoly(x,n)
    USE nrtype 
    implicit none 
    REAL(DP), DIMENSION(:), intent(in) :: x
    integeR(I4B), intent(in) :: n
    REAL(DP), Dimension(n) :: funpoly
    integer :: j
    REAL(DP) :: x1
    funpoly(1) = 1
    x1 = x(1) 
    do j = 2, n
       funpoly(j) = funpoly(j-1)*x1
    end do
    return
  end FUNCTION funpoly

  FUNCTION funre(x,n)
    USE nrtype 
    implicit none 
    REAL(DP), DIMENSION(:), intent(in) :: x
    integeR(I4B), intent(in) :: n
    REAL(DP), Dimension(n) :: funre
    REAL(DP) :: x1
    x1 = x(1) 
    funre(1) = 1.d0
    
    funre(2) = x1**(1.d0/4.d0)

    !funre(2) = x1**1.5d0
    funre(4)  = x1**0.75d0
    funre(3) = x1**0.5d0 
    return
  end FUNCTION funre

  
  FUNCTION funvol(x,n)
    USE nrtype 
    implicit none 
    REAL(DP), DIMENSION(:), intent(in) :: x
    integeR(I4B), intent(in) :: n
    REAL(DP), Dimension(n) :: funvol
    REAL(DP) :: x1
    x1 = x(1) 

    
!!$    funstokes(1) = x
!!$    funstokes(4) = x*x
!!$    funstokes(2) = x**1.5d0
!!$    funstokes(3) = x**(1.d0/3.d0)

    funvol(1) = (x1**1.d0)/((1.d0-x1)**3.d0)
    
    funvol(2) = x1**(2.d0)/((1.d0-x1)**3.d0)
    
    funvol(3)  = (x1**(1.5d0))/((1.d0-x1)**3.d0)
    funvol(4)  = (x1**(1.d0/3.d0))/((1.d0-x1)**3.d0)
    !/((1.d0-x1)**3.d0)
    !funvol(5)  = (x1**(0.57))!/((1.d0-x1)**3.d0)
    !funvol(4) = 1.d0
    !funvol(4) = 1.d0
    return
  end FUNCTION funvol

  FUNCTION funstokes(x,n)
    USE nrtype 
    implicit none 
    REAL(DP), DIMENSION(:), intent(in) :: x
    integeR(I4B), intent(in) :: n
    REAL(DP), Dimension(n) :: funstokes
    integer :: j, idim 
    REAL(DP) :: onebyeighteen, x1, x2, tmp

    idim = size(x)

    x1 = x(1) 
    x2 = 0.d0

    if(idim.eq.2) x2 = x(2) 
    tmp = 1.d0!/(1.d0-x1)**3.d0

    funstokes(4) = 0.d0!x2**2.d0*x1**1.5/(1-x1)**2.d0
    !funstokes(4) = x2*x1**2.0/(1-x1**1.5)**4d0 ! Works well
    funstokes(1) = x2*x1**6.0/(1-x1)**2.d0!**5.0!**3.d0!(1-dsqrt(x1))**3.0*x2!*x1*(1-x1)**2.d0 !Best
    !funstokes(2) = (x2**0.687)/(1-x1)**3.d0 ! Best
    funstokes(2) = x1**3.d0*x2
    
    funstokes(3) = (x1)/(1-x1)**3.d0 ! Works very well
    funstokes(5) = 0.0!x1**1.5*x2(1-x1)**4
    funstokes(6) =  x1**0.33/(1-x1)**4.d0
    
    !if(x2.gt.0.0d0) then 

       !Trials
       !funstokes(5) = x2***x1/(1-x1)**3.d0
    !else
!       funstokes(1) = 0.d0
!       funstokes(2) = 0.d0 
     !  WRITE(*,*)' RE = ', x2
     !  READ(*,*)

       !funstokes(3) = 0.d0
       !funstokes(5) = 0.d0 
    !end if
   !funstokes(6) = 0.d0 !x1**0.5d0 
    !funstokes(6) = (x1**(1.d0/3.d0))
    !WRITE(*,*)'x2 = ', x2, funstokes(3)

    !funstokes(6) = 0.d0
    !funstokes(5) = x2**0.4d0

    !else! if(x2.gt.5.d0.and.x2.le.10.d0) then 
       !.45-.5 good for fcc
       !funstokes(5) = 0.d0
       !funstokes(6) = ((1.d0-x1)**(0.d0))*(x2**0.35d0)/tmp
       !funstokes(4) = 0.d0

   ! else if(x2.gt.10.d0) then 
   !    funstokes(5) = 0.d0
   !    funstokes(4) = ((1.d0-x1)**(0.d0))*(x2**0.35d0)/tmp
   !    funstokes(6) = 0.d0
 
!    end if

    !funstokes(4) = (x1**(2.d0))

    !funstokes(4) = (x1**0.1d0)*((1.d0 - x1)**2.D0)*(x2**0.05d0)
    !funstokes(4) = ((x1)**2.D0)*(x2**0.05d0)


    !funstokes(4) = (x1**2.d0)
    !funstokes(5) = ((x1)**3.d0)*(exp(x2*0.001)) !x1*x1*x2!*log(x2)!*x1/(x2**1.d0)!/(x2**1.d0)

!!$
!!$    funstokes(1) = x1*x1
!!$    funstokes(2) = x1*x2
!!$    funstokes(3) = x2

    !funstokes(7) = 1.d0

    !funstokes(6) = ((1.d0-x1)**(3.d0))*(x2**0.8d0) !x1*x1*x2!*log(x2)!*x1/(x2**1.d0)!/(x2**1.d0)
    !funstokes(7) = ((0.7d0-x1)**(3.d0))*(x2**(0.1d0)) !x1*x1*x2!*log(x2)!*x1/(x2**1.d0)!/(x2**1.d0)
    !funstokes(7) = ((1.d0-x1)**(-4.d0))*(x2*x2) !x1*x1*x2!*log(x2)!*x1/(x2**1.d0)!/(x2**1.d0)
    !if(x1.eq.0.d0)WRITE(*,*) 'BASIS FUNCTIONS ARE', funstokes(:)
    !funstokes(6) = x2**(1.5d0)!*log(x2)!*x1/(x2**1.d0)!/(x2**1.d0)
    !funstokes(4) = 1.d0
    !funstokes(4) = x1*x1*log(x2)!*x1/(x2**1.d0)!/(x2**1.d0)

    !funstokes(2) = x2*x2*x1 !/((1.d0-x1)**3.d0)

    !funstokes(2) = x1*x1/(x2**1.d0)!/(x2**1.d0)
    !funstokes(3) = (x1)*x2 !/(1.d0-x1)

    ! funstokes(4) = exp(-x1)

    !    funstokes(3) = x1/((1.d0-x1)**4.d0)
    Return
  end FUNCTION funstokes

  FUNCTION funnuss(x,n)
    USE nrtype 
    implicit none 
    REAL(DP), DIMENSION(:), intent(in) :: x
    integeR(I4B), intent(in) :: n
    REAL(DP), Dimension(n) :: funnuss
    integer :: j, idim 
    REAL(DP) :: onebyeighteen, x1, x2, tmp

    idim = size(x)

    x1 = x(1) 
    x2 = 0.d0

    if(idim.eq.2) x2 = x(2) 
    tmp = 1.d0!/(1.d0-x1)**3.d0

    funnuss(4) = 0.d0!exp(x2*x1) 
    
      
    funnuss(1) = x2**2d0*x1**3d0*(1-x1)**2.0!x2**2*x1/(1-x1)!**3.d0 !Best
    !funnuss(2) = (x2**0.687)/(1-x1)**3.d0 ! Best
    funnuss(2) = 0.0!x2*x1*(1-x1)!**1.5!(x2**0.526)/(1-x1)**3.d0 
    
    
    if(x2.gt.0.0d0) then 
       funnuss(3) = ((x1)/(1-x1))!**2.d0!**3.d0 ! Works very well
       funnuss(5) = x1**0.33/(1-x1)**1.5 !0.d0!1.d0/(1-x1)**3.d0
       funnuss(6) = x2*x1**1.5/(1-x1)**1.5!**2.d0
       !Trials
       !funstokes(5) = x2***x1/(1-x1)**3.d0
    else
       !       funstokes(1) = 0.d0
       !       funstokes(2) = 0.d0 
       WRITE(*,*)' RE = ', x2
       READ(*,*)
       
       !funstokes(3) = 0.d0
       !funstokes(5) = 0.d0 
    end if
    Return
  end FUNCTION funnuss

end MODULE fitting_functions
