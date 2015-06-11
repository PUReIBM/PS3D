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


module functions
  USE precision 
  USE constants
  USE general_funcs 
  USE global_data
  !independent module 
  implicit none 
  
Contains 
  SUBROUTINE hermite_polynomial(n, moments,weights,abscissa)
    IMPLICIT NONE
    INTEGER, Intent(in) :: n
    REAL(prcn), Intent(in),Dimension(:) :: moments
    REAL(prcn), Intent(out),Dimension(:) :: weights,abscissa
    REAL(prcn),ALLOCATABLE,DIMENSION(:) :: nu, beta,a,b
    REAL(prcn),ALLOCATABLE,DIMENSION(:,:) :: sig,Z
    INTEGER :: i, j, imax,sizem,k,l
    REAL(prcn) :: bmin, temp
    sizem = SIZE(moments)
    Write(*,'(1A10,2x,1I6)')'SIZE OF MOMENTS ARRAY = ', sizem

    ALLOCATE(nu(sizem),beta(sizem-1),a(n),b(n),sig(2*n+1,2*n+1),Z(n,n))

    Do j = 1, 2*n
       nu(j) = zero
       if(j.lt.2*n)beta(j) = zero
    End Do
    
    Do j = 1, n
       a(j) = zero
       b(j) = zero
       Do i = 1, n
          z(j,i) = zero
       End Do
    End Do  
    
    Do j = 1, 2*n+1
       Do i = 1, 2*n+1
          sig(i,j) = zero
       End Do
    End Do
    Do j = 1, 2*n
       imax = FLOOR((j-1)/2.d0)
       Do i = 0, imax 
          temp = moments(j-2*i)*((-1)**i)*factorial(j-1)/(4.d0**i)
          temp = temp/(factorial(j-1-2*i) * factorial(i))
          nu(j) = nu(j) + temp
       End Do 
       if(j.lt.2*n)beta(j) = (j-1)/two
    End Do
    PRINT*,'NU = ', nu(:) 
    Do i = 2, 2*n+1
       sig(2,i) = nu(i-1)
    End Do 
    a(1) = nu(2)/nu(1)
    b(1) = zero
    Do k = 3,(n+1)
       Do l = k,(2*n-k+3)
          sig(k,l) = sig(k-1,l+1)-a(k-2)*sig(k-1,l)-b(k-2)*sig(k-2,l)+beta(l-1)*sig(k-1,l-1) 
       End Do
       a(k-1) = sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1) 
       b(k-1) = sig(k,k)/sig(k-1,k-1) 
    End Do
    PRINT*,'A = ', a(:) 
    PRINT*,'B = ', b(:) 
    bmin = MINVAL(b)
    IF(bmin.lt.zero)then
      Write(*,'(A25)')'MOMENTS IN WHEELER HERMITE ARE NOT REALIZABLE'
      STOP
    ENDIF
    
    Do i = 1, n-1
       z(i,i) = a(i)
       z(i,i+1) = DSQRT(b(i+1))
       z(i+1,i) = z(i,i+1)
    End Do
    z(n,n) = a(n)
!    do i = 1,n
!       PRINT*,'z = ', Z(i,1:n)
!    end do
    if(n.eq.1)then
      weights(1) = 1.0
      abscissa(1) = 0.0
    else if(n.eq.2)then
      weights(1) = 0.5
      weights(2) = 0.5
      abscissa(1) = -1.0
      abscissa(2) = 1.0
    else if(n.eq.3)then
      weights(1) = 0.16667
      weights(2) = 0.66667
      weights(3) = 0.16667
      abscissa(1) = -dsqrt(3.d0)
      abscissa(2) = 0.0
      abscissa(3) = dsqrt(3.d0)
    else if(n.eq.4)then
      weights(1) = 0.16667
      weights(2) = 0.66667
      weights(3) = 0.16667
      weights(4) = 0.1666
      abscissa(1) = -dsqrt(3.d0)
      abscissa(2) = 0.0
      abscissa(3) = dsqrt(3.d0)
      abscissa(4) = dsqrt(3.d0)
   endif
  END SUBROUTINE hermite_polynomial

!!$  SUBROUTINE form_poisson_inverse
!!$    IMPLICIT NONE
!!$    REAL(prcn) :: lower(mx1-1), upper(mx1-1), diagonal(mx1), top_right, bot_left
!!$    INTEGER :: i
!!$    
!!$    do i = 1, mx1
!!$       diagonal(i) = -two/dx2-w2(j,k)
!!$       if(i.lt.mx1)then
!!$          lower(i)=one/dx2
!!$          upper(i)=one/dx2
!!$       end if
!!$    end do
!!$    top_right = one/dx2
!!$    bot_left = one/dx2
!!$    CALL INVERT_CYCLIC_TRIDIAG(mx1,lower,diagonal,upper,top_right,bot_left,poisson_inverse)
!!$
!!$  END SUBROUTINE form_poisson_inverse

!!$  SUBROUTINE form_momentum_inverse(rks)
!!$    IMPLICIT NONE
!!$    INTEGER, Intent(in) :: rks
!!$    REAL(prcn) :: lower(mx1-1), upper(mx1-1), diagonal(mx1), top_right, bot_left
!!$    INTEGER :: i
!!$    
!!$    do i = 1, mx1
!!$       diagonal(i) = one + dt*vis*coef(rks,2)*(g(i)*two/dx2+w2(j,k))
!!$       if(i.lt.mx1)then
!!$          lower(i)= -g(i)*dt*vis*coef(rks,2)/dx2
!!$          upper(i)= -g(i)*dt*vis*coef(rks,2)/dx2
!!$       end if
!!$    end do
!!$    top_right = dt*vis*coef(rks,2)/dx2
!!$    bot_left = dt*vis*coef(rks,2)/dx2
!!$    
!!$    CALL INVERT_CYCLIC_TRIDIAG(mx1,lower,diagonal,upper,top_right,bot_left,momentum_inverse(1:mx1,1:mx1))
!!$    
!!$  END SUBROUTINE form_momentum_inverse
  
  
  
  FUNCTION factorial(N)  
    INTEGER :: factorial 

    Integer:: i,fact
    Integer, Intent(in):: N

    fact = 1
    
    Do i = 1,N
       fact = fact * i
    End Do
    factorial = fact
    RETURN
  END FUNCTION factorial
  subroutine maxzero(u1,u2,u3)
    Integer :: my2, mz
    complex(prcn), dimension(:,:), Intent(inout) ::  u1, u2, u3
    !complex(prcn) ::  u2(my2,mz)
    !complex(prcn) ::  u3(my2,mz)
    my2 = size(u1,1)
    mz = size(u1,2)

    u1(1:my2,mz/2+1) = czero 
    u1(my2,1:mz) = czero

    u2(1:my2,mz/2+1) = czero 
    u2(my2,1:mz) = czero

    u3(1:my2,mz/2+1) = czero 
    u3(my2,1:mz) = czero
  end subroutine maxzero
end module functions
