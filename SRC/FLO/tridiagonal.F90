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

module tridiagonal
  !	*****************************************************************
#include "ibm.h"

  USE precision 
  USE constants 
  USE global_data
  !USe workarrays, Only : gam 
  IMPLICIT NONE
contains 

  subroutine tridag(a,b,c,r,u,n)
    implicit none     
    integer, intent(in) ::  n
    
    real(prcn) :: gam(mx)
    
    integer :: i,j
    
    real(prcn), Dimension(:), Intent(in)  :: a(n),b(n),c(n),r(n)
    real(prcn), Dimension(:), Intent(out)  :: u(n)
    real(prcn) ::  bet
    do i=1,n
       u(i)=zero
    enddo
    gam = zero
    if(b(1).eq.zero)pause
    
    bet=b(1)
    u(1)=r(1)/bet
    
    do j=2,n
       
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       
       if(bet.eq.zero) then 
          bet=1.d-30 
          write(19,*) 'zero pivot'
          pause
       endif
       
       u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    enddo
!    do i = 1, n
!       Write(*,'4(2x,g17.8)') u(i)
!    end do
!    STOP
  end subroutine tridag


  subroutine tridag3(a,b,c,r1,r2,r3,u1,u2,u3,n)
    implicit none     
    integer, intent(in) ::  n
    
    real(prcn) :: gam(n)
    
    integer :: i,j
    
    real(prcn), Dimension(:), Intent(in)  :: a(n),b(n),c(n),r1(n),r2(n),r3(n)
    real(prcn), Dimension(:), Intent(out)  :: u1(n),u2(n),u3(n)
    real(prcn) ::  bet

    do i=1,n
       u1(i)=zero
       u2(i)=zero
       u3(i)=zero
    enddo
    gam = zero
    if(b(1).eq.zero)pause
    
    bet=b(1)

    u1(1)=r1(1)/bet
    u2(1)=r2(1)/bet
    u3(1)=r3(1)/bet

    do j=2,n
       
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       
       if(bet.eq.zero) then 
          bet=1.d-30 
          write(19,*) 'zero pivot'
          pause
       endif
       
       u1(j)=(r1(j)-a(j)*u1(j-1))/bet
       u2(j)=(r2(j)-a(j)*u2(j-1))/bet
       u3(j)=(r3(j)-a(j)*u3(j-1))/bet
    enddo
    do j=n-1,1,-1
       u1(j)=u1(j)-gam(j+1)*u1(j+1)
       u2(j)=u2(j)-gam(j+1)*u2(j+1)
       u3(j)=u3(j)-gam(j+1)*u3(j+1)
    enddo


!    if(gauss_u)then
!      do i = 1, n
!         Write(*,'9(2x,g17.8)') a(i), b(i), c(i),r2(i), u1(i),u2(i),u3(i)
!      end do
!    STOP
!   endif
   ! if(gauss_u)then
    !  do i = 1, n
    !     Write(*,'6(2x,g17.8)') r1(i), u1(i), r2(i), u2(i), r3(i), u3(i)
    !  end do
    !STOP
   !endif
  end subroutine tridag3

  subroutine mpi_tridag(a, b, c, r, u, n)
    IMPLICIT NONE
    INTEGER, Intent(in) :: n
    REAL(prcn), Intent(in) :: a(n), b(n), c(n), r(n)
    REAL(prcn), Intent(inout) :: u(n)
    INTEGER :: i
    REAL(prcn) :: at(n), bt(n), ct(n), rt(n,1), sol(n,1)
#if PARALLEL
    do i = 1, n
       at(i) = a(i)
       bt(i) = b(i)
       ct(i) = c(i)
       rt(i,1) = r(i)
    end do 
    if(I_AM_NODE_ZERO) at(1) = zero
    if(myid.eq.nproc-1) ct(n) = zero
    CALL parallel_partition(at(1:n),bt(1:n),ct(1:n),1,rt(1:n,1),sol(1:n,1),n)
    do i = 1, n
       u(i) = sol(i,1)
    end do

#endif
  end subroutine mpi_tridag

  subroutine mpi_tridag3(a, b, c, r1, r2, r3, u1, u2, u3, n)
    IMPLICIT NONE
    INTEGER, Intent(in) :: n
    REAL(prcn), Intent(in) :: a(n), b(n), c(n), r1(n), r2(n), r3(n)
    REAL(prcn), Intent(inout) :: u1(n), u2(n), u3(n)
    INTEGER :: i
    REAL(prcn) :: at(n), bt(n), ct(n), rt(n,3), sol(n,3)
#if PARALLEL
    do i = 1, n
       at(i) = a(i)
       bt(i) = b(i)
       ct(i) = c(i)
       rt(i,1) = r1(i)
       rt(i,2) = r2(i)
       rt(i,3) = r3(i)
    end do
 
    if(I_AM_NODE_ZERO) at(1) = zero
    if(myid.eq.nproc-1) ct(n) = zero

    CALL parallel_partition(at(1:n),bt(1:n),ct(1:n),3,rt(1:n,1:3),sol(1:n,1:3),n)
    do i = 1, n
       u1(i) = sol(i,1)
       u2(i) = sol(i,2)
       u3(i) = sol(i,3)
    end do

#endif
  end subroutine mpi_tridag3


  SUBROUTINE parallel_partition(a, b, c, nr, r, u, n)
    IMPLICIT NONE
    INTEGER, Intent(in) :: n, nr
    REAL(prcn), Intent(in) :: a(n), b(n), c(n)
    REAL(prcn), Intent(inout) :: r(n,nr)
    REAL(prcn), Intent(out) :: u(n,nr)
    REAL(prcn) :: alk, blk, clk, rlk(nr), auk, buk, cuk, ruk(nr), alphai, betai, ainf(2), binf(2), cinf(2), rinf(2,nr)
    REAL(prcn) :: inftri_a(2*nproc), inftri_b(2*nproc), inftri_c(2*nproc), inftri_r(2*nproc, nr), inftri_sol(2*nproc, nr)
    INTEGER :: i, ir
#if PARALLEL
    !PRINT*,' IN PARALLEL PARTITION'
    !Do i = 1, n
	!Write(*,'(6(2x,g17.8)') a(i),b(i),c(i),r(i,1:nr)
    !End do
    alk = a(2)
    blk = b(2)
    clk = c(2)
    rlk(1:nr) = r(2, 1:nr)
    ! Lower interface equation
    do i = 3, n
       alphai = a(i)/blk
       alk = -alphai*alk
       blk = b(i) - alphai*clk
       clk = c(i)
       do ir = 1, nr
          rlk(ir) = r(i,ir) - alphai*rlk(ir)
       end do
    end do
    ! Upper interface equation
    auk = a(n-1)
    buk = b(n-1)
    cuk = c(n-1)
    ruk(1:nr) = r(n-1, 1:nr)

    do i = n-2, 1, -1
       betai = c(i)/buk
       cuk = -betai*cuk
       buk = b(i) - betai*auk
       auk = a(i)
       do ir = 1, nr
          ruk(ir) = r(i,ir) - betai*ruk(ir)
       end do
    end do
    ainf(1) = auk
    ainf(2) = alk
    binf(1) = buk
    binf(2) = blk
    cinf(1) = cuk
    cinf(2) = clk
    rinf(1,1:nr) = ruk(1:nr)
    rinf(2,1:nr) = rlk(1:nr)
    !Write(*,'(A,I,6(2x,g17.8))')'INTERFACE EQUATIONS: ', myid, ainf(1), binf(1), cinf(1), rinf(1,1), rinf(1,2), rinf(1,3)
    
    CALL MPI_GATHER(ainf(1), 2, MPI_DOUBLE_PRECISION, inftri_a(1), 2, MPI_DOUBLE_PRECISION, NODE_ZERO, decomp_group, err_code)
    CALL MPI_GATHER(binf(1), 2, MPI_DOUBLE_PRECISION, inftri_b(1), 2, MPI_DOUBLE_PRECISION, NODE_ZERO, decomp_group, err_code)
    CALL MPI_GATHER(cinf(1), 2, MPI_DOUBLE_PRECISION, inftri_c(1), 2, MPI_DOUBLE_PRECISION, NODE_ZERO, decomp_group, err_code)
    do ir = 1, nr
       CALL MPI_GATHER(rinf(1,ir), 2, MPI_DOUBLE_PRECISION, inftri_r(1,ir), 2, MPI_DOUBLE_PRECISION, NODE_ZERO, decomp_group, err_code)
    end do
    if(I_AM_NODE_ZERO)then
     !do i = 1, 2*nproc
        !Write(*,'(6(2x,g17.8))')inftri_a(i), inftri_b(i), inftri_c(i), inftri_r(i,1), inftri_r(i,2), inftri_r(i,3)
     !end do 
       CALL multi_tridag(inftri_a, inftri_b, inftri_c, nr, inftri_r, inftri_sol, 2*nproc)
     !if(gauss_u)then
     !   do i = 1, 2*nproc
     !      Write(*,'(6(2x,g17.8)') inftri_sol(i,1:3)
     !   end do
       ! endif
    endif
   
     do ir = 1, nr
        CALL MPI_SCATTER(inftri_sol(1,ir), 2, MPI_DOUBLE_PRECISION, rinf(1,ir), 2, MPI_DOUBLE_PRECISION, NODE_ZERO, decomp_group, err_code)
     end do
     do ir = 1, nr
        u(1,ir) = rinf(1,ir)
        u(n,ir) = rinf(2,ir)
        r(2,ir) = r(2,ir) - a(2)*u(1,ir)
        r(n-1,ir) = r(n-1,ir) - c(n-1)*u(n,ir)
      end do
      CALL multi_tridag(a(2:n-1), b(2:n-1), c(2:n-1), nr, r(2:n-1,1:nr), u(2:n-1,1:nr), n-2)

    !do i = 1, n
    !   Write(*,'(A,I,6(2x,g17.8))')'SOLUTION: ', i, u(i,1:nr)
    !end do
    
#endif
  END SUBROUTINE parallel_partition

  subroutine multi_tridag(a,b,c,nr,r,u,n)
    implicit none     
    integer, intent(in) ::  n, nr
    
    real(prcn) :: gam(mx)
    
    integer :: i,j, ir
    
    real(prcn), Dimension(:), Intent(in)  :: a(n),b(n),c(n)
    real(prcn), Dimension(:,:), Intent(inout)  :: r(n,nr), u(n,nr)
    real(prcn) ::  bet
  
    do i=1,n
       u(i,1:nr)=zero
    enddo
    gam = zero
    if(b(1).eq.zero)pause
    
    bet=b(1)
    u(1,1:nr)=r(1,1:nr)/bet
    
    do j=2,n
       
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       
       if(bet.eq.zero) then 
          bet=1.d-30 
          write(19,*) 'zero pivot'
          pause
       endif
       do ir = 1, nr
          u(j,ir)=(r(j,ir)-a(j)*u(j-1,ir))/bet
       end do
    enddo
    do j=n-1,1,-1
       u(j,1:nr)=u(j,1:nr)-gam(j+1)*u(j+1,1:nr)
    enddo
!    do i = 1, n
!       Write(*,'4(2x,g17.8)') u(i)
!    end do
!    STOP
  end subroutine multi_tridag

  subroutine tridag_prefixsum(a, b, c, r, u, n)
    IMPLICIT NONE
    INTEGER, Intent(in) :: n
    REAL(prcn), Intent(in) :: a(n), b(n), c(n), r(n)
    REAL(prcn), Intent(inout) :: u(n)
    REAL(prcn) :: Bi(n+1,3,3), sendbuf(6), recvbuf(6), Brecv(3,3), xzero
    INTEGER :: i, j, sendproc, recvproc, ncomsteps, index, d, i1, j1
#if PARALLEL    
    do i = 1, 3
       do j = i,3
          if(i.eq.j)then
             Bi(1,i,j) = one
           !  Brecv(i,j) = one
          else
             Bi(1,i,j) = zero
             Bi(1,j,i) = Bi(1,i,j)
            ! Brecv(i,j) = zero
            ! Brecv(j,i) = Brecv(i,j)
          end if
       end do
    end do

    do i = 2, n+1
       Brecv(1,1) = -b(i-1)/c(i-1)
       Brecv(1,2) = -a(i-1)/c(i-1)
       Brecv(1,3) = r(i-1)/c(i-1)
       Brecv(2,1) = one
       Brecv(2,2) = zero
       Brecv(2,3) = zero

       Brecv(3,1) = zero
       Brecv(3,2) = zero
       Brecv(3,3) = one

       if(gauss_u)then
	 PRINT*,' I = ', i
         do d = 1, 3
            Write(*,'(A,3(2x,g12.5))')'Bij = ', Brecv(d,1:3)
         end do
       endif !STOP
       Do i1 = 1, 3
          Do j1 = 1, 3
             Bi(i,i1,j1) = Brecv(i1,1)*Bi(i-1,1,j1)+Brecv(i1,2)*Bi(i-1,2,j1)+Brecv(i1,3)*Bi(i-1,3,j1)
          End Do
       End Do
      if(gauss_u)then
      
       do d = 1, 3
          Write(*,'(A,3(2x,g12.5))')'Bij = ', Bi(i,d,1:3)
       end do
      endif
       !Bi(i,1:3,1:3) = MATMUL(Bi(i,1:3,1:3),Bi(i-1,1:3,1:3)) ! B_i = B_i B_{i-1}
       !if(i.eq.5)then       
       ! do d = 1, 3
       !   Write(*,'(A,3(2x,g12.5))')'Bij A = ', Bi(i,d,1:3)
       !end do
       !PARALLEL_FINISH()
       !STOP
       !endif
    end do
    ncomsteps = INT(log10(real(nproc,prcn))/log10(two))
    !PRINT*,'NCOMSTEPS = ', ncomsteps
    Do d = 0, ncomsteps
       do i = 1, 2
          do j = 1, 3
             index = (i-1)*3+j
             sendbuf(index) = Bi(n+1,i,j)
          end do
       end do
       sendproc = INT(2**d) + myid
       recvproc = myid - INT(2**d) 
       if(sendproc.gt.nproc-1) sendproc = MPI_PROC_NULL
       if(recvproc.lt.node_zero) recvproc = MPI_PROC_NULL
       !PRINT*, 'sendproc = ', sendproc
       !PRINT*, 'recvproc = ', recvproc
       RSENDRECV(sendbuf(1),6,sendproc,1,recvbuf(1),6,recvproc,1,decomp_group,status)
       if(.not.(recvproc.eq.MPI_PROC_NULL))then
         do i = 1, 2
            do j = 1, 3
               index = (i-1)*3+j
               Brecv(i,j) = recvbuf(index)
            end do
         end do
         do i = 1, n+1
            Bi(i,1:3,1:3) = MATMUL(Bi(i,1:3,1:3),Brecv(1:3,1:3)) ! B_i = B_i B_{i-1}
         end do
       end if
   End Do    
   If(myid.eq.nproc-1) then
     xzero = -Bi(n+1,1,3)/Bi(n+1,1,1)
   Endif
   BROADCAST_DOUBLE(xzero, 1, NPROC-1,comm_group)
   Do i = 1, n
      u(i) = Bi(i,1,1)*xzero + Bi(i,1,3)
   End Do 
#endif
  end subroutine tridag_prefixsum
end module tridiagonal
















