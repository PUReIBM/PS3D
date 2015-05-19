
MODULE randomno
  !Now it is a dependent module, RG 04/29/08
  USE global_data 
  USE precision            ! independent modules 
  USE constants
  Implicit None 
  Private
  
  real(prcn), Dimension(:), Allocatable, Save :: seedarray
  Integer, Save :: count, seedarraycount=400
  Public :: uni_dist, jn_dist, moment1, norm_dist
  Contains
    Subroutine norm_dist(frandn)
      Implicit None
      Integer :: N
      Real(prcn) , intent(out),dimension(:) :: frandn
      Real(prcn) , Dimension(:),allocatable :: frandtemp
      Integer :: sdvalue, i
      real(prcn) :: tempsdvalue
      Integer, Dimension(:), Allocatable :: rseed

      open(4,file = TRIM(RUN_NAME)//'_seedvalue', status = 'unknown')

      N=size(frandn)
      
      Allocate(frandtemp(N))
      call getnewseed(sdvalue)
      call ransed(-1,4,sdvalue,1)
      call rann2(frandtemp,1,1,1,N,1)
      frandn = frandtemp
      deallocate(frandtemp)
    End Subroutine norm_dist

    Subroutine uni_dist(frandu)
      Implicit None
      integer:: N
      real(prcn),intent(out),dimension(:) :: frandu
      real(prcn), Dimension(:), Allocatable :: frandtemp
      integer :: sdvalue,i
      real(prcn) :: tempsdvalue
      Integer, Dimension(:), Allocatable :: rseed
      
      open(4,file = TRIM(RUN_NAME)//'_seedvalue', status = 'unknown')
      N=size(frandu)
      !Print*,'in uni dist'
      Allocate(frandtemp(N))
      call getnewseed(sdvalue)
      
      !Print*,'in uni dist: after get new seed'
      call ransed(-1,4,sdvalue,1)
      call ranu2(frandtemp,1,1,1,N,1)
      frandu = frandtemp
      deallocate(frandtemp)
    End Subroutine uni_dist
    
    Subroutine jn_dist(frandn,np,ndim,umf0,rsf)
      Implicit None
      Integer, Intent(in) :: np,ndim
      real(prcn),intent(out),dimension(:,:) :: frandn
      Real(prcn), Intent(in),dimension(:,:) ::  rsf
      Real(prcn), Intent(in),dimension(:) :: umf0
      real(prcn), dimension(:,:), Allocatable :: frandntemp
      Real(prcn), dimension(:,:), Allocatable ::  rsftemp
      Real(prcn), dimension(:), Allocatable :: umftemp
      Integer :: sdvalue
      Allocate(frandntemp(size(frandn,1),size(frandn,2)))
      Allocate(rsftemp(size(rsf,1),size(rsf,2)))
      Allocate(umftemp(size(umf0,1)))
      
      rsftemp = rsf
      umftemp = umf0
      open(4,file = TRIM(RUN_NAME)//'_seedvalue', status = 'unknown')
      !print*,'np = ',np,' ndim = ',ndim,' rsf=',rsf,' umf0=',umf0
      
      call getnewseed(sdvalue)
      call ransed(-1,4,sdvalue,1)
      
      call ranjn(frandntemp,1,np,1,ndim,np,umftemp,rsftemp,ndim)
      frandn = frandntemp
      Deallocate(frandntemp, rsftemp, umftemp)
    End Subroutine jn_dist
    
     Subroutine getnewseed(sdvalue)
      Implicit None
      Integer, Intent(out) :: sdvalue
      Integer :: tempseed,testsdvalue, trycount
      Logical :: change_seed
      change_seed = .false.
      If(.not.Allocated(seedarray)) Allocate(seedarray(seedarraycount))
      count = count + 1
      trycount = 0 
      !print*,'in randomno.... count = ',count
16    if (count.gt.seedarraycount.or.count.lt.1) count = 1
      trycount = trycount + 1
      IF(TRYCOUNT.GT.20) PRINT*,'FAILING IN GETNEWSEED: ALREADY TRIED 20 COMPLETE SEED ARRAY COUNTS'
      If(count.eq.1) then
         change_seed = .true.
         Open(1001,file= TRIM(RUN_NAME)//'_seed.d')
         Read(1001,39) tempseed
         call ransed(-1,4,tempseed,1)
         call ranu2(seedarray,1,1,1,seedarraycount,1)
         close(1001)
      end If
39    Format((I8))
18    sdvalue = seedarray(count)*100000000
      testsdvalue=sdvalue/10000000
      !print*,'in randomno.... count = ',testsdvalue, sdvalue
            
      If (testsdvalue.eq.0) then
         count = count + 1
         If(count.lt.seedarraycount) Then 
            go to 18
         Else 
            !count = 1
            go to 16
         End If
      End If
      !rewind(1)
      !print*,'New value of seed = ', sdvalue
      If(change_seed) then 
         open(1001, file = TRIM(RUN_NAME)//'_seed.d')
         write(1001,39) sdvalue
         !write(1,39) 12345678!sdvalue
         close(1001)
      endif
    End Subroutine getnewseed
!!$c___________________________________________________________
!!$c
!!$c     Subroutine to calculate the moments of a
!!$c     sample : cf Numerical Recipes.
!!$c     modified for UNEQUAL weights.
!!$c     modified so it calculates only the moments needed.
!!$c     mode = 1 mean only
!!$c     mode = 2 mean and variance
!!$c     mode = 3 mean, variance and skewness
!!$c     mode = 4 mean, variance, skewness and kurtosis.
!!$c
!!$c___________________________________________________________

!    call moment1(4, n, n, wt, u, 1e-06, ave,adev,sdev,svar,skew&
   subroutine moment1( mode, nd, n, w, data,ave, adev, sdev&
        &,var, skew, curt )
     implicit none
     integer , intent(in) :: nd,n, mode
     real(prcn), intent(in), Dimension(:) ::  data,w
     !Real(prcn), Intent(in) :: vlmt
     Real(prcn), Intent(out) :: ave, adev, sdev, var, skew,&
          & curt
     real(prcn) ::  dsum,dwt,dfluc,dvar,dskew,dcurt, dflucw, vlmt
     Integer :: j
     if(n.le.1)pause 'ERROR in moment: N must be atleast 2'
!

     vlmt=1e-025
     dsum=0.d0
     adev=0.
     var=0.
     skew=0.
     curt=0.

     dwt = 0.d0
     dfluc = 0.d0
     dflucw = 0.d0
     dsum= 0.d0
     dvar=0.d0
     dskew=0.d0
     dcurt=0.d0

     do j = 1,n
       dsum = dsum + (w(j)*data(j))
     enddo
     ave=(dsum)

     if ( mode .eq. 1 ) then
       return
     elseif( mode .eq. 2 )then
       do j = 1,n
         dwt = (w(j))
         dfluc=(data(j)-ave)
         dflucw=dfluc*dwt
         dsum = dsum + abs(dflucw)
         dflucw = dflucw*dfluc
         dvar=dvar+ dflucw
       enddo
       adev=dsum
       var = dvar
       sdev=sqrt(var)
     elseif( mode .eq. 3 )then
       do j = 1,n
         dwt = (w(j))
         dfluc=(data(j)-ave)
         dflucw=dfluc*dwt
         dsum = dsum + abs(dflucw)
         dflucw = dflucw*dfluc
         dvar=dvar+ dflucw
         dflucw = dflucw*dfluc
         dskew=dskew+dflucw
       enddo
       adev=(dsum)
       var = (dvar)
       skew= (dskew)
       sdev=sqrt(var)
       if(var.gt.vlmt) then
         skew=skew/(sdev**3)
       else
         write(*,*)'WARNING in moment: Variance less than vlmt'
       endif
     elseif( mode .eq. 4 )then
       do j = 1,n
         dwt = (w(j))
         dfluc=(data(j)-ave)
         dflucw=dfluc*dwt
         dsum = dsum + abs(dflucw)
         dflucw = dflucw*dfluc
         dvar=dvar+ dflucw
         dflucw = dflucw*dfluc
         dskew=dskew+dflucw
         dflucw = dflucw*dfluc
         dcurt=dcurt+dflucw
       enddo
       adev=(dsum)
       var = (dvar)
       skew= (dskew)
       curt = (dcurt)
       sdev=sqrt(var)
       if(var.gt.vlmt) then
         skew=skew/(sdev**3)
         curt=curt/(var**2)!- 3.
       endif
     endif
     return
   end subroutine moment1


  End MODULE randomno


SUBROUTINE RANLUX(RVEC,LENV)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!          
!       references: 
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++ 
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	DIMENSION RVEC(LENV)
	DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
	PARAMETER (MAXLEV=4, LXDFLT=3)
	DIMENSION NDSKIP(0:MAXLEV)
	DIMENSION NEXT(24)
	PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
	PARAMETER (ITWO24=2**24, ICONS=2147483563)
	SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
	SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
	INTEGER LUXLEV
	LOGICAL NOTYET
	DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
	DATA I24,J24,CARRY/24,10,0./

!Luxury Level   0     1     2   *3*    4
	DATA NDSKIP/0,   24,   73,  199,  365 /
!Corresp. to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
	IF (NOTYET) THEN
		NOTYET = .FALSE.
		JSEED = JSDFLT  
		INSEED = JSEED
		WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
		LUXLEV = LXDFLT
		NSKIP = NDSKIP(LUXLEV)
		LP = NSKIP + 24
		IN24 = 0
		KOUNT = 0
		MKOUNT = 0
		WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ', LUXLEV,'      p =',LP
		TWOM24 = 1.
		DO I= 1, 24
			TWOM24 = TWOM24 * 0.5
			K = JSEED/53668
			JSEED = 40014*(JSEED-K*53668) -K*12211
			IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
			ISEEDS(I) = MOD(JSEED,ITWO24)
		ENDDO
		TWOM12 = TWOM24 * 4096.
		DO I= 1,24
			SEEDS(I) = REAL(ISEEDS(I))*TWOM24
			NEXT(I) = I-1
		ENDDO
		NEXT(1) = 24
		I24 = 24
		J24 = 10
		CARRY = 0.
		IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
	ENDIF

!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
	DO IVEC= 1, LENV
		UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
		IF (UNI .LT. 0.)  THEN
			UNI = UNI + 1.0
			CARRY = TWOM24
		ELSE
			CARRY = 0.
		ENDIF
		SEEDS(I24) = UNI
		I24 = NEXT(I24)
		J24 = NEXT(J24)
		RVEC(IVEC) = UNI
!  small numbers (with less than 12 "significant" bits) are "padded".
		IF (UNI .LT. TWOM12)  THEN
			RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
!        and zero is forbidden in case someone takes a logarithm
			IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
		ENDIF
!        Skipping to luxury.  As proposed by Martin Luscher.
		IN24 = IN24 + 1
		IF (IN24 .EQ. 24)  THEN
			IN24 = 0
			KOUNT = KOUNT + NSKIP
			DO ISK= 1, NSKIP
				UNI = SEEDS(J24) - SEEDS(I24) - CARRY
				IF (UNI .LT. 0.)  THEN
					UNI = UNI + 1.0
					CARRY = TWOM24
				ELSE
					CARRY = 0.
				ENDIF
				SEEDS(I24) = UNI
				I24 = NEXT(I24)
				J24 = NEXT(J24)
			ENDDO
		ENDIF
	ENDDO
	KOUNT = KOUNT + LENV
	IF (KOUNT .GE. IGIGA)  THEN
		MKOUNT = MKOUNT + 1
		KOUNT = KOUNT - IGIGA
	ENDIF
	RETURN

!           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
	TWOM24 = 1.
	DO 195 I= 1, 24
		NEXT(I) = I-1
	195 TWOM24 = TWOM24 * 0.5
	NEXT(1) = 24
	TWOM12 = TWOM24 * 4096.
	WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
	WRITE(6,'(5X,5I12)') ISDEXT
	DO I= 1, 24
		SEEDS(I) = REAL(ISDEXT(I))*TWOM24
	ENDDO
	CARRY = 0.
	IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
	ISD = IABS(ISDEXT(25))
	I24 = MOD(ISD,100)
	ISD = ISD/100
	J24 = MOD(ISD,100)
	ISD = ISD/100
	IN24 = MOD(ISD,100)
	ISD = ISD/100
	LUXLEV = ISD
	IF (LUXLEV .LE. MAXLEV) THEN
		NSKIP = NDSKIP(LUXLEV)
		WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',LUXLEV
	ELSE  IF (LUXLEV .GE. 24) THEN
		NSKIP = LUXLEV - 24
		WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
	ELSE
		NSKIP = NDSKIP(MAXLEV)
		WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
		LUXLEV = MAXLEV
	ENDIF
	INSEED = -1
	RETURN
!
!                    Entry to ouput seeds as integers
	ENTRY RLUXUT(ISDEXT)
	DO I= 1, 24
		ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
	ENDDO
	ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
	IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
	RETURN
!
!                   Entry to output the "convenient" restart point
	ENTRY RLUXAT(LOUT,INOUT,K1,K2)
	LOUT = LUXLEV
	INOUT = INSEED
	K1 = KOUNT
	K2 = MKOUNT
	RETURN
!
!                   Entry to initialize from one or three integers
	ENTRY RLUXGO(LUX,INS,K1,K2)
	IF (LUX .LT. 0) THEN
		LUXLEV = LXDFLT
	ELSE IF (LUX .LE. MAXLEV) THEN
		LUXLEV = LUX
	ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
		LUXLEV = MAXLEV
		WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
	ELSE
		LUXLEV = LUX
		DO ILX= 0, MAXLEV
			IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
		ENDDO
	ENDIF
	IF (LUXLEV .LE. MAXLEV)  THEN
		NSKIP = NDSKIP(LUXLEV)
		WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',LUXLEV,'     P=', NSKIP+24
	ELSE
		NSKIP = LUXLEV - 24
		WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
	ENDIF
	IN24 = 0
	IF (INS .LT. 0)  WRITE (6,'(A)') ' Illegal initialization by RLUXGO, negative input seed'
	IF (INS .GT. 0)  THEN
		JSEED = INS
		WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',JSEED, K1,K2
	ELSE
		JSEED = JSDFLT
		WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
	ENDIF
	INSEED = JSEED
	NOTYET = .FALSE.
	TWOM24 = 1.
	DO I= 1, 24
		TWOM24 = TWOM24 * 0.5
		K = JSEED/53668
		JSEED = 40014*(JSEED-K*53668) -K*12211
		IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
		ISEEDS(I) = MOD(JSEED,ITWO24)
	ENDDO
	TWOM12 = TWOM24 * 4096.
	DO I= 1,24
		SEEDS(I) = REAL(ISEEDS(I))*TWOM24
		NEXT(I) = I-1
	ENDDO
	NEXT(1) = 24
	I24 = 24
	J24 = 10
	CARRY = 0.
	IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
	KOUNT = K1
	MKOUNT = K2
	IF (K1+K2 .NE. 0)  THEN
		DO IOUTER= 1, K2+1
			INNER = IGIGA
			IF (IOUTER .EQ. K2+1)  INNER = K1
			DO ISK= 1, INNER
				UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
				IF (UNI .LT. 0.)  THEN
					UNI = UNI + 1.0
					CARRY = TWOM24
				ELSE
					CARRY = 0.
				ENDIF
				SEEDS(I24) = UNI
				I24 = NEXT(I24)
				J24 = NEXT(J24)
			ENDDO
		ENDDO
!         Get the right value of IN24 by direct calculation
		IN24 = MOD(KOUNT, NSKIP+24)
		IF (MKOUNT .GT. 0)  THEN
			IZIP = MOD(IGIGA, NSKIP+24)
			IZIP2 = MKOUNT*IZIP + IN24
			IN24 = MOD(IZIP2, NSKIP+24)
		ENDIF
!       Now IN24 had better be between zero and 23 inclusive
		IF (IN24 .GT. 23) THEN
			WRITE (6,'(A/A,3I11,A,I5)') '  Error in RESTARTING with RLUXGO:','  The values', INS, &
			&    K1, K2, ' cannot occur at luxury level', LUXLEV
			IN24 = 0
		ENDIF
	ENDIF
	RETURN
END
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
