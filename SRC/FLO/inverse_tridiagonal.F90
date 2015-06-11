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


MODULE inverse_tridiagonal
#include "ibm.h"
  USE precision
  USE constants
  IMPLICIT NONE
CONTAINS
  
SUBROUTINE INVERT_TRIDIAG(n,lower,diag,upper,Ainv)
  IMPLICIT NONE
  INTEGER, Intent(in) :: n
  REAL*8, Intent(out) :: Ainv(1:n,1:n)
  REAL*8,Intent(in) :: lower(1:n-1), diag(1:n),upper(1:n-1)
  REAL*8 :: theta(0:n), phi(n+1)
  INTEGER :: i,j, ii,jj
  PRINT*,' INVERTING THE TRIDIAGONAL MATRIX'
  
  ! Initial conditions on theta
  theta(0) = 1
  theta(1) = diag(1)
  ! Initial conditions on phi
  phi(n) = diag(n)
  phi(n+1) = 1

  do i = 2,n
     theta(i) = diag(i)*theta(i-1) - upper(i-1)*lower(i-1)*theta(i-2)
  end do

  PRINT*,'SMALL NUMBER IS : ', SMALL_NUMBER

  if(ABS(theta(n)).le.SMALL_NUMBER)then
     PRINT*,'SINGULAR MATRIX. INVERSE DOES NOT EXIST.'
  end if

  do i = n-1,1,-1
     phi(i) = diag(i)*phi(i+1) - lower(i)*upper(i)*phi(i+2)
  end do
  
  do i = 1, n
     do j = 1, n
        Ainv(i,j) = 1.D0
        if(i.le.j)then
           do ii = i,j-1
              Ainv(i,j) = Ainv(i,j) * upper(ii)
           end do
           Ainv(i,j) = Ainv(i,j)*theta(i-1)*phi(j+1)/theta(n)
        else if(i.gt.j)then 
           do jj = j,i-1
              Ainv(i,j) = Ainv(i,j)*lower(jj)
           end do
           Ainv(i,j) = Ainv(i,j)*theta(j-1)*phi(i+1)/theta(n)
        end if
        Ainv(i,j) = (-1)**(i+j)*Ainv(i,j)
     end do
  end do
END SUBROUTINE INVERT_TRIDIAG

SUBROUTINE INVERT_CYCLIC_TRIDIAG(n,lower,diag,upper,top_right,bot_left,Ainv)
  IMPLICIT NONE
  INTEGER, Intent(in) :: n
  REAL(prcn), Intent(out) :: Ainv(1:n,1:n)
  REAL(prcn),Intent(in) :: lower(1:n-1), diag(1:n),upper(1:n-1),top_right, bot_left
  REAL(prcn) :: diag_tri(1:n), Ainv_tri(n,n)
  REAL(prcn) :: alpha, beta, u1, un, v1, vn, vtainvu, opvtainvu, Apij
  INTEGER :: i, j

  alpha = one
  beta = one
  do i = 1, n
     diag_tri(i) = diag(i)
  end do
  u1 = alpha
  un = beta
  
  v1 = bot_left/un
  vn = top_right/u1

  diag_tri(1) = diag_tri(1) - u1*v1
  diag_tri(n) = diag_tri(n) - un*vn

  CALL invert_tridiag(N,lower(1:N-1),diag_tri(1:N),upper(1:N-1),Ainv_tri(1:N,1:N))
  vtainvu = zero
  vtainvu = v1*(Ainv_tri(1,1)*u1 + Ainv_tri(1,n)*un)
  vtainvu = vtainvu + vn*(Ainv_tri(n,1)*u1 + Ainv_tri(n,n)*un)
  opvtainvu = one + vtainvu
  IF(ABS(opvtainvu).le.SMALL_NUMBER)then
    Write(*,'(A50)')'SINGULAR MATRIX. INVERSE DOES NOT EXIST. STOPPING'
    PARALLEL_FINISH()
    STOP
  ENDIF
  do i = 1, n 
     do j = 1, n
        Apij = (Ainv_tri(i,1)*u1*v1+Ainv_tri(i,n)*un*v1)*Ainv_tri(1,j)
        Apij = Apij + (Ainv_tri(i,1)*u1*vn+Ainv_tri(i,n)*un*vn)*Ainv_tri(n,j)
        Ainv(i,j) = Ainv_tri(i,j) - Apij/opvtainvu
     end do
  end do
END SUBROUTINE INVERT_CYCLIC_TRIDIAG
END MODULE inverse_tridiagonal
