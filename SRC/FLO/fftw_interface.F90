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


MODULE fftw_interface
#include "ibm.h"
  !Independent Module ....
  USE precision             ! independent modules
  USE constants
  USE global_data  !, ONLY : my, mz, my2
  IMPLICIT NONE
  SAVE
  PUBLIC

CONTAINS
#if !FFTW3
  SUBROUTINE ff2rc(a,b)
    IMPLICIT NONE
    !REAL to COMPLEX fft

    REAL(prcn), INTENT(in), DImension(:,:)  :: a
    COMPLEX(prcn), INTENT(out), DImension(:,:)  :: b

    INTEGER*8 :: plan
    INTEGER :: arr(2)
    INTEGER :: i,j,k, mylo, mzlo, my2lo !Local values of my, mz, my2
    mylo = size(a,1)
    mzlo = size(a,2)
    my2lo = size(b,1)
!!$    Print*, 'mylo = ', mylo,mzlo,my2lo
!!$
!!$     print*,'integer size = ',bit_size(plan)
!!$     read(*,*)
    arr(1)=mylo
    arr(2)=mzlo

    CALL rfftwnd_f77_create_plan(plan,2,arr,-1,0)
    CALL rfftwnd_f77_one_real_to_complex(plan,a,b)

    DO j=1,mzlo,1
       DO i=1,my2lo,1
          b(i,j)=b(i,j)/(float(mylo*mzlo))
       ENDDO
    ENDDO


    CALL rfftwnd_f77_destroy_plan(plan)

    RETURN
  END SUBROUTINE ff2rc

  SUBROUTINE ff2cr(a,b)
    !COMPLEX to REAL fft
    IMPLICIT NONE
    SAVE
    COMPLEX(prcn), INTENT(in), DIMENSION(:,:)  :: a

    REAL(prcn), INTENT(out) , DIMENSION(:,:) :: b

    COMPLEX(prcn), DIMENSION(:,:), ALLOCATABLE, SAVE  :: a1
    INTEGER*8 ::  plan
    INTEGER ::  arr(2)

    Integer, save :: count=0
    INTEGER :: i,j,k, mylo, mzlo, my2lo !Local values of my, mz, my2

    mylo = size(b,1)
    mzlo = size(b,2)
    my2lo = size(a,1)

    count  = count + 1
    If(count.eq.1) allocate(a1(my2lo,mzlo))
!!$    Print*, 'mylo = ', mylo,mzlo,my2lo
!!$
!!$    print*,'integer size = ',bit_size(plan)

    arr(1)=mylo
    arr(2)=mzlo

    DO j=1,mzlo,1
       DO i=1,my2lo,1
          a1(i,j)=a(i,j)
       ENDDO
    ENDDO

    CALL rfftwnd_f77_create_plan(plan,2,arr,1,0)
    CALL rfftwnd_f77_one_complex_to_real(plan,a1,b)

!!$    DO j=1,mz,1
!!$       DO i=1,my,1
!!$          b(i,j)=b(i,j)/(float(my*mz))
!!$       ENDDO
!!$    ENDDO

    CALL rfftwnd_f77_destroy_plan(plan)

    if(count.gt.10000000) then !To make sure over flow does occur for
       ! count variable
       count = 0
       deallocate(a1)
    end if
  END SUBROUTINE ff2cr
!!$
!!$subroutine ff2rc1(a,b)
!!$
!!$
!!$    real(prcn) ::  a(my,mz)
!!$    complex(prcn) :: b(my2,mz)
!!$
!!$    integer*8 :: plan
!!$    integer ::  arr(2)
!!$    integer :: i,j,k
!!$    !plan =
!!$    print*,'integer size = ',bit_size(plan)
!!$    read(*,*)
!!$    arr(1)=my
!!$    arr(2)=mz
!!$    !Print*,'my, mz =', my,mz,' a(2:2) = ', a(1:2,2:4), b
!!$    call rfftwnd_f77_create_plan(plan,2,arr,-1,0)
!!$    call rfftwnd_f77_one_real_to_complex(plan,a,b)
!!$
!!$    do j=1,mz,1
!!$       do i=1,my2,1
!!$          b(i,j)=b(i,j)/(float(my*mz))
!!$       enddo
!!$    enddo
!!$
!!$
!!$    call rfftwnd_f77_destroy_plan(plan)
!!$
!!$    return
!!$  end subroutine ff2rc1
!!$
      SUBROUTINE ff3rc(freal,fhat)

      COMPLEX(prcn), INTENT(out),DIMENSION(:,:,:) :: fhat
      REAL(prcn), INTENT(in),DIMENSION(:,:,:) ::  freal
      REAL(prcn),ALLOCATABLE,DIMENSION(:,:,:), SAVE:: frealtemp
      INTEGER, SAVE :: count=0
      INTEGER*8 ::  plan
      INTEGER :: arr(3)
      INTEGER ::  i,j,k

     !dimension arr(ndim)
      arr(1)=size(freal,1)
      arr(2)=size(freal,2)
      arr(3)=size(freal,3)
      count = count+1
      If(count.eq.1) allocate(frealtemp(arr(1),arr(2),arr(3)))

      DO k=1,arr(3)
         DO j=1,arr(2)
            DO i=1,arr(1)
               frealtemp(i,j,k)=freal(i,j,k)
            ENDDO
         ENDDO
      ENDDO

     CALL rfftwnd_f77_create_plan(plan,3,arr,-1,0)
     CALL rfftwnd_f77_one_real_to_complex(plan,frealtemp,fhat)

      DO k=1,arr(3)
         DO j=1,arr(2)
            DO i=1,arr(1)/2+1
               fhat(i,j,k)=fhat(i,j,k)/(DBLE(arr(1)*arr(2)*arr(3)))
            ENDDO
         ENDDO
      ENDDO

      CALL rfftwnd_f77_destroy_plan(plan)
      if(count.gt.10000000) then !To make sure over flow does occur for
                                ! count variable
         count = 0
         deallocate(frealtemp)
      end if
      RETURN
      END SUBROUTINE ff3rc

!///////////////////////////////////////////////////////////////////
!       complex to real fft for 3D array

	SUBROUTINE ff3cr(fhat,freal)

		COMPLEX(prcn),INTENT(in),DIMENSION(:,:,:) ::  fhat
		REAL(prcn),INTENT(in),DIMENSION(:,:,:)  ::  freal
		COMPLEX(prcn),ALLOCATABLE,DIMENSION(:,:,:), SAVE :: fhattemp
		INTEGER, SAVE :: count=0
		INTEGER*8 ::  plan
		INTEGER :: arr(3)
		INTEGER :: i,j,k

		arr(1)=size(freal,1)
		arr(2)=size(freal,2)
		arr(3)=size(freal,3)

		!^^^^^^^ Mohammad 10-27-2009 ^^^^^^^
		count = count+1
		If(count.eq.1) allocate(fhattemp(arr(1)/2+1,arr(2),arr(3)))
		!-----------------------------------
		DO k=1,arr(3)
			DO j=1,arr(2)
				DO i=1,arr(1)/2+1
					fhattemp(i,j,k)=fhat(i,j,k)
				ENDDO
			ENDDO
		ENDDO

		CALL rfftwnd_f77_create_plan(plan,3,arr,1,0)
		CALL rfftwnd_f77_one_complex_to_real(plan,fhattemp,freal)
		CALL rfftwnd_f77_destroy_plan(plan)
		if(count.gt.10000000) then !To make sure over flow does occur for
			! count variable
			count = 0
			deallocate(fhattemp)
		end if
		RETURN
      END SUBROUTINE ff3cr

#else 

  SUBROUTINE ff2rc(a,b)
    IMPLICIT NONE 
    !REAL to COMPLEX fft
    
    REAL(prcn), INTENT(in), DImension(:,:)  :: a
    COMPLEX(prcn), INTENT(out), DImension(:,:)  :: b
    INTEGER*8 :: planr2c_2d
    REAL*8 ,ALLOCATABLE,DIMENSION(:,:),SAVE :: frealtemp
      Integer, save :: count=0
    INTEGER :: i,j,k, mylo, mzlo, my2lo !Local values of my, mz, my2

    mylo = size(a,1)
    mzlo = size(a,2)
    my2lo = size(b,1)
!    Print*, 'mylo = ', mylo,mzlo,my2lo
!!$    
!!$     print*,'integer size = ',bit_size(plan)
!$$$      read(*,*)
    count  = count + 1 
    If(count.eq.1)  ALLOCATE(frealtemp(mylo,mzlo))
    
      do j=1,mzlo
         do i=1,mylo
            frealtemp(i,j)=a(i,j)
         enddo
      enddo
      call dfftw_plan_dft_r2c_2d(planr2c_2d,mylo,mzlo,frealtemp,b,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(planr2c_2d, frealtemp, b)
      call dfftw_destroy_plan(planr2c_2d)
    DO j=1,mzlo,1
       DO i=1,my2lo,1
          b(i,j)=b(i,j)/(float(mylo*mzlo))               
       ENDDO
    ENDDO
      if(count.gt.10000000) then !To make sure over flow does occur for
                                ! count variable 
         count = 0
         DEALLOCATE(frealtemp)
      end if
    
    RETURN
  END SUBROUTINE ff2rc

  SUBROUTINE ff2cr(a,b)
    !COMPLEX to REAL fft
    IMPLICIT NONE 
    SAVE 
    COMPLEX(prcn), INTENT(in), DIMENSION(:,:)  :: a
    
    REAL(prcn), INTENT(out) , DIMENSION(:,:) :: b
    
    COMPLEX(prcn), DIMENSION(:,:), ALLOCATABLE, SAVE  :: fhattemp

    INTEGER*8 :: planc2r_2d

    Integer, save :: count 
    INTEGER :: i,j,k, mylo, mzlo, my2lo !Local values of my, mz, my2
    
    mylo = size(b,1)
    mzlo = size(b,2)
    my2lo = size(a,1)
    
    count  = count + 1 
    If(count.eq.1) allocate(fhattemp(my2lo,mzlo))
!!$    Print*, 'mylo = ', mylo,mzlo,my2lo
!!$    
!!$    print*,'integer size = ',bit_size(plan)

    
    DO j=1,mzlo,1
       DO i=1,my2lo,1
          fhattemp(i,j)=a(i,j)
       ENDDO
    ENDDO
      call dfftw_plan_dft_c2r_2d(planc2r_2d,mylo,mzlo,fhattemp,b,FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(planc2r_2d, fhattemp, b)
      call dfftw_destroy_plan(planc2r_2d)
      if(count.gt.10000000) then !To make sure over flow does occur for
                                ! count variable 
         count = 0
         DEALLOCATE(fhattemp)
      end if
      END SUBROUTINE ff2cr
#endif 
END MODULE fftw_interface
