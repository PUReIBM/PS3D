!-------
! Independent module
!-------
! Some general purpose functions and routines
!-------
! Author: Chidambaram Narayanan
!         Nuclear Engineering Laboratory
!         ETH Zurich
! Author: RAHUL GARG
!         IOWA STATE UNIVERSITY
!-------
MODULE general_funcs
#include "../FLO/ibm.h"
  use global_data
  USE precision
  USE constants

  Implicit none
  Private
  Public:: mag, locate, hunt, blankline, separator, getnewunit &
          ,average_value, root_mean_square, array_dot_product, screen_separator, cross_product3, openfile, gener_filename, cross_product, eigen_val_vec, instant_file_opener

  !-------
  ! Generic interface
  !-------
  INTERFACE average_value
    Module procedure average_oned_array
    Module procedure average_twod_array
    Module procedure average_threed_array
  END INTERFACE

  INTERFACE root_mean_square
    Module procedure rms_oned_array
    Module procedure rms_twod_array
    Module procedure rms_threed_array
  END INTERFACE

  INTERFACE array_dot_product
    Module procedure array_dot_product_1d
    Module procedure array_dot_product_2d
    Module procedure array_dot_product_3d
  END INTERFACE

!-------
 Contains
!-------

	SUBROUTINE GENER_FILENAME(FILENAME, FILENAMEIN)
		use global_data
		IMPLICIT NONE

		Character(LEN=*),Intent(out) :: filename
		Character(LEN=*),Intent(in) :: FILENAMEIN

		Character*80 :: FILENAMELOC
		Integer :: node, strlen

		if(I_AM_NODE_ZERO)then
			FILENAME = TRIM(FILENAMEIN)
			FILENAMELOC = ""
		else
			FILENAME = ""
		end if

#if PARALLEL
		if(I_AM_NODE_ZERO)then
			do node=1,nproc-1
				FILENAMELOC = TRIM(FILENAMEIN)

				if (.not.communicator_done) then
					SEND_STRING(filenameloc,strlen,node,0,1,comm_group)
				else
					SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
				endif
			end do
		else  
			if (.not.communicator_done) then
				RECV_STRING(filename,strlen,node_zero,0,1,comm_group,status)
			else
				RECV_STRING(filename,strlen,node_zero,0,1,decomp_group,status)
			endif
		endif
#endif
	END SUBROUTINE GENER_FILENAME


	subroutine instant_file_opener(filename_in,unitnum,formated)
		implicit none
		character*100, intent(in) :: filename_in
		integer, intent(in) :: unitnum
		logical, intent(in) :: formated
		character*100 filename, type
		logical :: filexist

		if (formated) then
			type="formatted"
		else
			type="unformatted"
		endif

		if (from_post) then
			filename=trim(filename_in)//"_post.dat"
			open (unitnum, file=trim(filename), status="replace", action="write", form=trim(type))
			if (formated) write (unitnum,*) "zone"
		else
			filename=trim(filename_in)//".dat"
			if (irestart==0.and.first_pass) then
				open (unitnum, file=trim(filename), status="replace", action="write", form=trim(type))
				if (formated) write (unitnum,*) "zone"
			else
				inquire (file=trim(filename), exist=filexist)
				if (.not.filexist) then
					write (*,*) 'THE FILE "'//trim(filename)//'" DOES NOT EXIST. OPENNING A NEW ONE'
					open (unitnum, file=trim(filename), status="replace", action="write", form=trim(type))
					if (formated) write (unitnum,*) "zone"
				else
					open (unitnum, file=trim(filename), status="old", action="write", form=trim(type), position="append")
				endif
			endif
		endif
	end subroutine instant_file_opener



  !----------
  !----------
  ! A . B where A and B are arrays of the same dimensions
  !----------
  !----------
  FUNCTION array_dot_product_3d(A,B)  !3D Arrays
    Real(prcn):: array_dot_product_3d

    Integer:: i, j, k, s1, s2, s3
    Real(prcn), Dimension(:,:,:), Intent(in):: A, B

    s1 = SIZE(A,1)
    s2 = SIZE(A,2)
    s3 = SIZE(A,3)

    array_dot_product_3d = zero
    Do 10 i = 1,s1
    Do 10 j = 1,s2
    Do 10 k = 1,s3
      array_dot_product_3d = array_dot_product_3d   &
                            + A(i,j,k)*B(i,j,k)
    10 Continue
  END FUNCTION array_dot_product_3d

  FUNCTION array_dot_product_2d(A,B)  !2D Arrays
    Real(prcn):: array_dot_product_2d

    Integer:: i, j, s1, s2
    Real(prcn), Dimension(:,:), Intent(in):: A, B

    s1 = SIZE(A,1)
    s2 = SIZE(A,2)

    array_dot_product_2d = zero
    Do 10 i = 1,s1
    Do 10 j = 1,s2
      array_dot_product_2d = array_dot_product_2d   &
                            + A(i,j)*B(i,j)
    10 Continue
  END FUNCTION array_dot_product_2d

  FUNCTION array_dot_product_1d(A,B)  !1D Arrays
    Real(prcn):: array_dot_product_1d

    Integer:: i, s1
    Real(prcn), Dimension(:), Intent(in):: A, B

    s1 = SIZE(A,1)

    array_dot_product_1d = zero
    Do 10 i = 1,s1
      array_dot_product_1d = array_dot_product_1d   &
                            + A(i)*B(i)
    10 Continue
  END FUNCTION array_dot_product_1d
  
  !----------
  !----------
  ! Calculate the rms of a one-dimensional array
  !----------
  !----------
  FUNCTION rms_oned_array(onedarray,knownavg)
    Real(prcn):: rms_oned_array
    Real(prcn), Intent(in), Optional:: knownavg
    Real(prcn), Dimension(:):: onedarray

    Integer:: n, npts
    Real(prcn):: average, rms

    npts = SIZE(onedarray)

    If (PRESENT(knownavg)) then
      average = knownavg
    Else
      average = average_value(onedarray)
    Endif

    rms = 0.0
    Do n = 1,npts
      rms = rms + (onedarray(n) - average)**2
    Enddo

    rms_oned_array = SQRT(rms/REAL(npts))
  END FUNCTION rms_oned_array

  !----------
  !----------
  ! Calculate the rms of a two-dimensional array
  !----------
  !----------
  FUNCTION rms_twod_array(twodarray,knownavg)
    Real(prcn):: rms_twod_array
    Real(prcn), Intent(in), Optional:: knownavg
    Real(prcn), Dimension(:,:):: twodarray

    Integer:: m, n, mpts, npts
    Real(prcn):: average, rms

    mpts = SIZE(twodarray,1)
    npts = SIZE(twodarray,2)

    If (PRESENT(knownavg)) then
      average = knownavg
    Else
      average = average_value(twodarray)
    Endif

    rms = 0.0
    Do 10 m = 1,mpts
    Do 10 n = 1,npts
      rms = rms + (twodarray(m,n) - average)**2
    10 Continue

    rms_twod_array = SQRT(rms/(REAL(mpts)*REAL(npts)))
  END FUNCTION rms_twod_array

  !----------
  !----------
  ! Calculate the rms of a three-dimensional array
  !----------
  !----------
  FUNCTION rms_threed_array(threedarray,knownavg)
    Real(prcn):: rms_threed_array
    Real(prcn), Intent(in), Optional:: knownavg
    Real(prcn), Dimension(:,:,:):: threedarray

    Integer:: el, m, n, elpts, mpts, npts
    Real(prcn):: average, rms

    elpts = SIZE(threedarray,1)
    mpts  = SIZE(threedarray,2)
    npts  = SIZE(threedarray,3)

    If (PRESENT(knownavg)) then
      average = knownavg
    Else
      average = average_value(threedarray)
    Endif

    rms = 0.0
    Do 10 el = 1,elpts
    Do 10 m  = 1,mpts
    Do 10 n  = 1,npts
      rms = rms + (threedarray(el,m,n) - average)**2
    10 Continue

    rms_threed_array = SQRT(rms  &
                     /(REAL(elpts)*REAL(mpts)*REAL(npts)))
  END FUNCTION rms_threed_array

  !----------
  !----------
  ! Calculate the average value of a one dimensional array
  !----------
  !----------
  FUNCTION average_oned_array(onedarray)
    Real(prcn):: average_oned_array
    Real(prcn), Dimension(:):: onedarray

    Integer:: n, npts
    Real(prcn):: avg_value

    npts = SIZE(onedarray)

    avg_value = 0.0
    Do n = 1,npts
      avg_value = avg_value + onedarray(n)
    Enddo
    average_oned_array = avg_value/REAL(npts)
  END FUNCTION average_oned_array

  !----------
  !----------
  ! Calculate the average value of a two dimensional array
  !----------
  !----------
  FUNCTION average_twod_array(twodarray)
    Real(prcn):: average_twod_array
    Real(prcn), Dimension(:,:):: twodarray

    Integer:: m, n, mpts, npts
    Real(prcn):: avg_value

    mpts = SIZE(twodarray,1)
    npts = SIZE(twodarray,2)

    avg_value = 0.0
    Do 10 m = 1,mpts
    Do 10 n = 1,npts
      avg_value = avg_value + twodarray(m,n)
    10 Continue

    average_twod_array = avg_value/(REAL(mpts)*REAL(npts))
  END FUNCTION average_twod_array

  !----------
  !----------
  ! Calculate the average value of a three dimensional array
  !----------
  !----------
  FUNCTION average_threed_array(threedarray)
    Real(prcn):: average_threed_array
    Real(prcn), Dimension(:,:,:):: threedarray

    Integer:: el, m, n, elpts, mpts, npts
    Real(prcn):: avg_value

    elpts = SIZE(threedarray,1)
    mpts  = SIZE(threedarray,2)
    npts  = SIZE(threedarray,3)

    avg_value = 0.0
    Do 10 el = 1,elpts
    Do 10 m  = 1,mpts
    Do 10 n  = 1,npts
      avg_value = avg_value + threedarray(el,m,n)
    10 Continue

    average_threed_array = avg_value/(REAL(elpts)*REAL(mpts)*REAL(npts))
  END FUNCTION average_threed_array

  !----------
  !----------
  ! Pick an unused unit number
  !----------
  !----------
  SUBROUTINE openfile(unit, filename, formfile, statusfile)
    Implicit none 
    Character(LEN=*), optional :: filename, formfile, statusfile
    INTEGER, INTENT(out) :: unit 
    
    LOGICAL :: filexist, isopen
   
    IF(.not.present(statusfile)) then 
       statusfile = 'unknown'
       WRITE(*,'(A,A)')'NO STATUS SPECIFIED FOR FILE: ', filename

       WRITE(*,'(A)') ' THEREFORE OPENING A FILE WITH STATUS = unknown'
    end IF
    
    IF(.not.present(formfile)) formfile = 'formatted'

! WRITE(*,'(A,/,A,/A,/,A)')'OPENING FILENAME WITH FORM AND STATUS AS  ', filename, formfile, statusfile
    OPEN(unit, file=TRIM(filename), form=TRIM(formfile), status=TRIM(statusfile), err=1000)

    RETURN
1000  CONTINUE
    WRITE(*,'(A,A)') 'COULD NOT OPEN FILE: ', filename
    
    INQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)

    WRITE(*,*) 'DOES THIS FILE EXIST ?', filexist

    WRITE(*,*) 'IS IT ALREADY OPEN ?', isopen
    
    RETURN 
  end SUBROUTINE openfile


  FUNCTION getnewunit(minval,maxval)
    Integer:: getnewunit
    Integer, Intent(in):: minval, maxval
   
    Integer:: iu, ios
    Logical:: isokay, isused

    Do iu = minval, maxval
      INQUIRE(UNIT=iu, EXIST=isokay, OPENED=isused &
             ,IOSTAT=ios)
      If (isokay.and.(.not.isused).and.(ios.eq.0)) then
        getnewunit = iu
        return
      Endif
    Enddo
    getnewunit = -1
  END FUNCTION getnewunit

  !----------
  !----------
  ! To put a line of dashes in a formatted file
  !----------
  !----------
  SUBROUTINE CROSS_PRODUCT3 (AA, XX,YY) 
    IMPLICIT NONE
    ! 
    REAL(prcn):: AA(3), XX(3), YY(3) 
    ! 
    AA(1) = XX(2)*YY(3) - XX(3)*YY(2) 
    AA(2) = XX(3)*YY(1) - XX(1)*YY(3) 
    AA(3) = XX(1)*YY(2) - XX(2)*YY(1)

    RETURN  
  END SUBROUTINE CROSS_PRODUCT3


  SUBROUTINE separator(unitno,length,symbol)
    Integer, Intent(in):: unitno, length
    Character(LEN=1):: symbol

    Integer:: i

    Do i = 1,length
      write(unitno,'(A1)',ADVANCE='NO')symbol
    Enddo
    write(unitno,'()',ADVANCE='YES')
  END SUBROUTINE separator

  SUBROUTINE screen_separator(length,symbol)
    Integer, Intent(in):: length
    Character(LEN=1):: symbol

    Integer:: i

    Do i = 1,length
      write(*,'(A1)',ADVANCE='NO')symbol
    Enddo
    write(*,'()',ADVANCE='YES')
  END SUBROUTINE screen_separator


  !----------
  !----------
  ! To put a blank line in a formatted file
  !----------
  !----------
  SUBROUTINE blankline(unitno)
    Integer, Intent(in):: unitno

    write(unitno,*)' '
  END SUBROUTINE blankline

  !----------
  !----------
  ! To calculate the magnitude of a vector
  !----------
  !----------
  FUNCTION mag(u)
    Real(prcn):: mag
    Real(prcn), Dimension(:):: u

    Integer:: i, length

    length = SIZE(u)

    mag = 0.0
    Do i = 1,length
      mag = mag + u(i)**2.0
    Enddo
    mag = sqrt(mag)
  END FUNCTION mag

  !----------
  !----------
  ! Bisection method to search an array
  ! Array 'x' can be in either ascending or descending order
  ! Note: New logical operator .eqv. being used
  ! If xval is out of range then locate is either 0 or n
  ! Lifted from Numerical Recipies in Fortran 90
  !----------
  !----------
  FUNCTION locate(x,xval)
    Integer:: locate
    Real(prcn), Dimension(:), Intent(in):: x
    Real(prcn), Intent(in):: xval

    Integer:: n, jl, jm, ju
    Logical:: ascend

    n = size(x)
    ascend = (x(n) >= x(1))

    jl = 0
    ju = n+1

    If (xval == x(1)) then
      locate = 1
      return
    Elseif (xval == x(n)) then
      locate = n-1
      return
    Endif

    Do
      if ((ju - jl) <= 1) exit
      jm = (ju + jl)/2
      if (ascend .eqv. (xval >= x(jm))) then
        jl = jm
      else
        ju = jm
      endif
    Enddo
    locate = jl
  END FUNCTION locate

  !----------
  !----------
  ! Hunt and then bisect
  ! Starts search from a given guess value
  ! Lifted from Numerical Recipies in Fortran 90
  !----------
  !----------
  FUNCTION hunt(x,xval,iguess)
    Integer:: hunt
    Real(prcn), Dimension(:), Intent(in):: x
    Real(prcn), Intent(in):: xval
    Integer:: iguess

    Integer:: n, inc, jhi, jlo, jm
    Logical:: ascend
    !print*,'in hunt n=',n
    n = size(x)
    jlo = iguess
    ascend = (x(n) >= x(1))

    IF (jlo <= 0 .or. jlo > n) then
      jlo = 0
      jhi = n+1
    ELSE
      inc = 1
      If (xval >= x(jlo) .eqv. ascend) then
        Do
          jhi = jlo + inc
          if (jhi > n) then
            jhi = n+1
            exit
          else
            if (xval < x(jhi) .eqv. ascend) exit
            jlo = jhi
            inc = inc + inc
          endif
        Enddo
      Else
        jhi = jlo
        Do 
          jlo = jhi - inc
          if (jlo < 1) then
            jlo = 0
            exit
          else
            if (xval >= x(jlo) .eqv. ascend) exit
            jhi = jlo
            inc = inc + inc
          endif
        Enddo
      Endif
    ENDIF

    Do
      If (jhi-jlo <= 1) then
        if (xval == x(n)) jlo = n-1
        if (xval == x(1)) jlo = 1
        exit
      Else
        jm = (jhi + jlo)/2
        if (xval >= x(jm) .eqv. ascend) then
          jlo = jm
        else
          jhi = jm
        endif
      Endif
    Enddo
    hunt = jlo
  END FUNCTION hunt

	subroutine eigen_val_vec(a, n, w, vectors)
		implicit none
		integer, intent(in) :: n
		real(prcn), intent(in) :: a(n,n)
		real(prcn), intent(out) :: w(n), vectors(n,n)

		integer :: nselect, lda, ldz, info, lwork, liwork, il, iu, m, lwmax
		integer, allocatable :: isuppz(:), iwork(:)

		real(prcn) :: abstol, vl, vu
		real(prcn), allocatable :: work(:)
		real(prcn), allocatable :: tensor(:,:) 

!		EXTERNAL DSYEVR

		lwmax = 26*n
		ABSTOL = -1
		lda = n
		ldz = n

		allocate(isuppz(2*n), iwork(lwmax), work(lwmax), tensor(n,n))
		tensor = a

		IL = 1
		IU = NSELECT

		LWORK = -1
		LIWORK = -1

!		CALL DSYEVR('V', 'A', 'U', N, tensor, LDA, VL, VU, IL, IU, ABSTOL, M, W, vectors, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

		LWORK = MIN(LWMAX, INT(WORK(1)))
		LIWORK = MIN(LWMAX, IWORK(1))

		!Solve eigenproblem.
!		CALL DSYEVR('V', 'A', 'U', N, tensor, LDA, VL, VU, IL, IU, ABSTOL, M, W, vectors, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

		!Check for convergence.
		IF( INFO.GT.0 ) THEN
			WRITE(*,*)'The algorithm failed to compute eigenvalues.'
			STOP
		END IF
		deallocate(isuppz, iwork, work, tensor)
	end subroutine eigen_val_vec

	subroutine cross_product(vec1, vec2, vec3)
		implicit none
		real(prcn), intent(in) :: vec1(ndim), vec2(ndim)
		real(prcn), intent(out) :: vec3(ndim)

		vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
		vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
		vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
	end subroutine cross_product

END MODULE general_funcs

!-------
! List of functions
!-------
! o getnewunit - returns unused unit number
! o separator - puts a line of symbols in a file
! o blank - introduces a blank line in a file
! o mag - magnitude of a vector
! o locate - locate a particle 
! o hunt - locate a particle given a guess
!-------
