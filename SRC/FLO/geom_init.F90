MODULE geom_init
#include "ibm.h"
  USE precision            ! independent modules 
  USE constants
  USE global_data 
  IMPLICIT NONE   
CONTAINS
	SUBROUTINE quadgener(r,dr,f1,f2)
!!$
!!$c     ********************************************************************
!!$c     generate one quadrant of a spherical surface.
!!$c     dy,dz are assumed equal.
!!$c     ********************************************************************
!!$      
		INTEGER ::  imax,jmax
		INTEGER ::  i,j
		INTEGER ::  nbndloc,nrprloc
		REAL(prcn) ::  x,y,z,dtheta,theta,phi,dphi
		REAL(prcn) ::  xi,xo,yi,yo,zi,zo, xo2, yo2, zo2 
		REAL(prcn), INTENT(in) ::  r,f1,f2, dr
		REAL(prcn) ::  r2,ri, ro, ro2, rado
		REAL(prcn),DIMENSION(:), ALLOCATABLE :: nxb,nyb,nzb
		REAL(prcn),DIMENSION(:), ALLOCATABLE :: nxi,nyi,nzi
		REAL(prcn),DIMENSION(:), ALLOCATABLE :: nxo,nyo,nzo
		REAL(prcn),DIMENSION(:), ALLOCATABLE :: nxo2,nyo2,nzo2
!!$c     ********************************************************************
!!$c     generate surface points
!!$c     ********************************************************************
		nbnd=0
		nrpr=0
		theta=0.
		phi=0.
		nbnd=0

		!********************************************************************
		! generate flow reversal pairs
		!     goto 100

		ri=r-dr
		ro=r+dr
		ro2 = r+two*dr

		theta=0.
		phi=0.
		imax=NINT(twopi*r/four*f2)

		!imax = 30
		dtheta=twopi/imax/four
		nbndloc = 0
		nrprloc = 0

		DO i=1,imax+1
			jmax=NINT(twopi*r*SIN(theta)/4*f2)
			!jmax = 30
			dphi=twopi/jmax/4
			!Print*,'dphi  = ', dphi, r*SIN(theta)*SIN(dphi)
			phi=0.
			DO j=1,jmax+1
				rado = xo2**2. + yo2**2. + zo2**2.
				r2=x*x+y*y+z*z

				nbndloc=nbndloc+1
				phi=phi+dphi
				nrprloc=nrprloc+1
			ENDDO
			!read(*,*)
			theta=theta+dtheta
		ENDDO
		!PRINT*,'nbndloc 1= ', nbndloc
		ALLOCATE(nxb(nbndloc), nyb(nbndloc), nzb(nbndloc))
		!ALLOCATE(nxi(nbndloc), nyi(nbndloc), nzi(nbndloc))
		!ALLOCATE(nxo(nbndloc), nyo(nbndloc), nzo(nbndloc))
		!ALLOCATE(nxo2(nbndloc), nyo2(nbndloc), nzo2(nbndloc))

		ri=r-dr
		ro=r+dr
		ro2 = r+two*dr
    
		theta=0.
		phi=0.
		imax=NINT(twopi*r/four*f2)

		nbndloc = 0
		nrprloc = 0

		!Print*,'dtheta  = ', dtheta,imax
		DO i=1,imax+1
			jmax=NINT(twopi*r*SIN(theta)/4*f2)
			!jmax = 30
			dphi=twopi/jmax/4
			!Print*,'dphi  = ', dphi, r*SIN(theta)*SIN(dphi)

			phi=0.
			DO j=1,jmax+1
				!^^^^^^^^^^^ Mohammad 08-18-2010 ^^^^^^^^^^^^^^^^^^^^^^^^^		
				! changing the spherical coordinates direction to have the 
				! axis of the sphere aligned with the flow direction
				z=r*COS(theta)/r
				IF(ABS(z).LT.1.0d-5) z=0.
				y=r*SIN(theta)*COS(phi)/r
				IF(ABS(y).LT.1.0d-5) y=0.
				x=r*SIN(theta)*SIN(phi)/r
				IF(ABS(x).LT.1.0d-5) x=0.
				r2=x*x+y*y+z*z

				zi=ri*COS(theta)
				IF(ABS(zi).LT.1.0d-5) zi=0.
				yi=ri*SIN(theta)*COS(phi)
				IF(ABS(yi).LT.1.0d-5) yi=0.
				xi=ri*SIN(theta)*SIN(phi)        
				IF(ABS(xi).LT.1.0d-5) xi=0.
				zo=ro*COS(theta)
				IF(ABS(zo).LT.1.0d-5) zo=0.
				yo=ro*SIN(theta)*COS(phi)
				IF(ABS(yo).LT.1.0d-5) yo=0.
				xo=ro*SIN(theta)*SIN(phi)        
				IF(ABS(xo).LT.1.0d-5) xo=0.

				zo2=ro2*COS(theta)
				IF(ABS(zo2).LT.1.0d-5) zo2=0.
				yo2=ro2*SIN(theta)*COS(phi)
				IF(ABS(yo2).LT.1.0d-5) yo2=0.
				xo2=ro2*SIN(theta)*SIN(phi)        
				IF(ABS(xo2).LT.1.0d-5) xo2=0.

!		x=r*COS(theta)/r
!		IF(ABS(x).LT.1.0d-5) x=0.
!		y=r*SIN(theta)*COS(phi)/r
!		IF(ABS(y).LT.1.0d-5) y=0.
!		z=r*SIN(theta)*SIN(phi)/r
!		IF(ABS(z).LT.1.0d-5) z=0.
!		r2=x*x+y*y+z*z
!
!		xi=ri*COS(theta)
!		IF(ABS(xi).LT.1.0d-5) xi=0.
!		yi=ri*SIN(theta)*COS(phi)
!		IF(ABS(yi).LT.1.0d-5) yi=0.
!		zi=ri*SIN(theta)*SIN(phi)        
!		IF(ABS(zi).LT.1.0d-5) zi=0.
!
!		xo=ro*COS(theta)
!		IF(ABS(xo).LT.1.0d-5) xo=0.
!		yo=ro*SIN(theta)*COS(phi)
!		IF(ABS(yo).LT.1.0d-5) yo=0.
!		zo=ro*SIN(theta)*SIN(phi)        
!		IF(ABS(zo).LT.1.0d-5) zo=0.
!
!		xo2=ro2*COS(theta)
!		IF(ABS(xo2).LT.1.0d-5) xo2=0.
!		yo2=ro2*SIN(theta)*COS(phi)
!		IF(ABS(yo2).LT.1.0d-5) yo2=0.
!		zo2=ro2*SIN(theta)*SIN(phi)        
!		IF(ABS(zo2).LT.1.0d-5) zo2=0.

				rado = xo2**2. + yo2**2. + zo2**2.
				r2=x*x+y*y+z*z

				nbndloc=nbndloc+1
				phi=phi+dphi
				nrprloc=nrprloc+1
				nzb(nbndloc) = z 
				nyb(nbndloc) = y 
				nxb(nbndloc) = x 
!!$          nzi(nbndloc) = zi 
!!$          nyi(nbndloc) = yi 
!!$          nxi(nbndloc) = xi 
!!$          nzo(nbndloc) = zo 
!!$          nyo(nbndloc) = yo 
!!$          nxo(nbndloc) = xo 
!!$
!!$          nzo2(nbndloc) = zo2 
!!$          nyo2(nbndloc) = yo2 
!!$          nxo2(nbndloc) = xo2

			ENDDO
			!read(*,*)
			theta=theta+dtheta
		ENDDO
    
		!PRINT*,'nbndloc 2= ', nbndloc
		CALL generbnd(nbnd, nbndloc, nxb(1:nbndloc), nyb(1:nbndloc), nzb(1:nbndloc))
    
		DEALLOCATE(nxb, nyb, nzb)
	END SUBROUTINE quadgener

  SUBROUTINE generbnd(icount, npt, nx, ny, nz)

    IMPLICIT NONE

    REAL(prcn),DIMENSION(:), INTENT(IN) :: nx,ny,nz
    INTEGER, INTENT(IN) :: npt
    INTEGER :: i, nbnd_est, j
    INTEGER, INTENT(out) :: icount

    REAL(prcn) ::  cdloc(3)
    nbnd_est = 8*npt
    !if(I_AM_NODE_ZERO) OPEN(unit=11,file=TRIM(RUN_NAME)//'_bnd.dat',status='unknown')

    !PRINT*,'NBND_EST = ', nbnd_est, ndim, icount
    ALLOCATE(xs(ndim, nbnd_est),cd(ndim, nbnd_est))
    !ALLOCATE(nx(npt), ny(npt), nz(npt), da(npt))
    !READ(11,*) (nx(i),ny(i),nz(i),da(i),i=1,npt)

    icount = 0

    DO i=1,npt
       CALL dcgen(nx(i),ny(i),nz(i),cdloc)
       
       icount = icount + 1
       xs(1,icount) = nx(i)
       xs(2,icount) = ny(i)
       xs(3,icount) = nz(i)
       cd(1:3, icount) = cdloc(1:3)
       !WRITE(17,111) nx(i),ny(i),nz(i),da(i)
       !WRITE(18,112) cd(1),cd(2),cd(3)
    ENDDO

    DO i=1,npt
       IF(nx(i).NE.0) THEN
          CALL dcgen(-nx(i),ny(i),nz(i),cdloc)

          icount = icount + 1
          xs(1,icount) = -nx(i)
          xs(2,icount) = ny(i)
          xs(3,icount) = nz(i)
          cd(1:3, icount) = cdloc(1:3)
          !WRITE(17,111) -nx(i),ny(i),nz(i),da(i)
          !WRITE(18,112) cd(1),cd(2),cd(3)
       ENDIF
    ENDDO

    DO i=1,npt
       IF(ny(i).NE.0) THEN
          CALL dcgen(nx(i),-ny(i),nz(i),cdloc)

          icount = icount + 1     
          xs(1,icount) = nx(i)
          xs(2,icount) = -ny(i)
          xs(3,icount) = nz(i)
          cd(1:3, icount) = cdloc(1:3)
       ENDIF
    ENDDO

    DO i=1,npt
       IF(nz(i).NE.0) THEN
          CALL dcgen(nx(i),ny(i),-nz(i),cdloc)

          icount = icount + 1
          xs(1,icount) = nx(i)
          xs(2,icount) = ny(i)
          xs(3,icount) = -nz(i)
          cd(1:3, icount) = cdloc(1:3)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(ny(i).NE.0)) THEN
          CALL dcgen(-nx(i),-ny(i),nz(i),cdloc)

          icount = icount + 1
          xs(1,icount) = -nx(i)
          xs(2,icount) = -ny(i)
          xs(3,icount) = nz(i)
          cd(1:3, icount) = cdloc(1:3)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(nz(i).NE.0)) THEN
          CALL dcgen(-nx(i),ny(i),-nz(i),cdloc)

          icount = icount + 1
          xs(1,icount) = -nx(i)
          xs(2,icount) = ny(i)
          xs(3,icount) = -nz(i)
          cd(1:3, icount) = cdloc(1:3)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          CALL dcgen(nx(i),-ny(i),-nz(i),cdloc)
          icount = icount + 1
          xs(1,icount) = nx(i)
          xs(2,icount) = -ny(i)
          xs(3,icount) = -nz(i)
          cd(1:3, icount) = cdloc(1:3)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          CALL dcgen(-nx(i),-ny(i),-nz(i),cdloc)
          icount = icount + 1
          xs(1,icount) = -nx(i)
          xs(2,icount) = -ny(i)
          xs(3,icount) = -nz(i)
          cd(1:3, icount) = cdloc(1:3)

       ENDIF
    ENDDO
    !if(I_AM_NODE_ZERO)then
    !   do i = 1, icount 
    !      write(11,112) (xs(j,i), j = 1, 3)
    !   end do
    !   close(11, status="keep")
    !end if
112 FORMAT(5(E20.10,1x))
    
    !PRINT*,'TOTAL NUMNER ...ICOUTN = ', ICOUNT, size(xs,2)
  END SUBROUTINE generbnd


  SUBROUTINE generrpr(icount)

    IMPLICIT NONE

    INTEGER ::  npt, i,j
    INTEGER, INTENT(out) ::  icount
    REAL(prcn),DIMENSION(:), ALLOCATABLE :: nx,ny,nz, c1
    REAL(prcn) :: jacobian(3,3)
    REAL(prcn) ::  jacobinv(3,3)

    OPEN(unit=11,file='sphrpi',status='old')
    OPEN(unit=12,file='sphrpo',status='old')
    OPEN(unit=14,file='sphrpo2',status='old')
    OPEN(unit=13,file='rpr',status='unknown')
    OPEN(unit=16,file='rprout',status='old')

    WRITE(*,*) 'Reading number of reversal points from rprout'
    READ(16,*) npt

    ALLOCATE(nx(npt), ny(npt), nz(npt), c1(npt))

    !First do it for the innner reversal points 
    READ(11,*) (nx(i),ny(i),nz(i),c1(i),i=1,npt)

    icount = 0

    DO i=1,npt
       WRITE(13,112) nx(i),ny(i),nz(i),c1(i)
       icount = icount + 1
    ENDDO

    DO i=1,npt
       IF(nx(i).NE.0) THEN
          WRITE(13,112) -nx(i),ny(i),nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO

    DO i=1,npt
       IF(ny(i).NE.0) THEN
          WRITE(13,112) nx(i),-ny(i),nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO

    DO i=1,npt
       IF(nz(i).NE.0) THEN
          WRITE(13,112) nx(i),ny(i),-nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO

    DO i=1,npt

       IF((nx(i).NE.0).AND.(ny(i).NE.0)) THEN

          WRITE(13,112) -nx(i),-ny(i),nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,112) -nx(i),ny(i),-nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO

    DO i=1,npt
       IF((ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,112) nx(i),-ny(i),-nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO

    DO i=1,npt

       IF((nx(i).NE.0).AND.(ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,112) -nx(i),-ny(i),-nz(i),c1(i)
          icount = icount + 1
       ENDIF
    ENDDO
    !Now do it for the first layer of outer points 
    READ(12,*) (nx(i),ny(i),nz(i),i=1,npt)

    DO i=1,npt
       WRITE(13,111) nx(i),ny(i),nz(i)
    ENDDO

    DO i=1,npt

       IF(nx(i).NE.0) THEN
          WRITE(13,111) -nx(i),ny(i),nz(i)
       ENDIF
    ENDDO

    DO i=1,npt

       IF(ny(i).NE.0) THEN

          WRITE(13,111) nx(i),-ny(i),nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF(nz(i).NE.0) THEN
          WRITE(13,111) nx(i),ny(i),-nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(ny(i).NE.0)) THEN
          WRITE(13,111) -nx(i),-ny(i),nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,111) -nx(i),ny(i),-nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,111) nx(i),-ny(i),-nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,111) -nx(i),-ny(i),-nz(i)
       ENDIF
    ENDDO


    !Finally do it for the third layer of  outer points RG 11/15/05
    READ(14,*) (nx(i),ny(i),nz(i),i=1,npt)

    DO i=1,npt
       WRITE(13,111) nx(i),ny(i),nz(i)
    ENDDO

    DO i=1,npt

       IF(nx(i).NE.0) THEN
          WRITE(13,111) -nx(i),ny(i),nz(i)
       ENDIF
    ENDDO

    DO i=1,npt

       IF(ny(i).NE.0) THEN

          WRITE(13,111) nx(i),-ny(i),nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF(nz(i).NE.0) THEN
          WRITE(13,111) nx(i),ny(i),-nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(ny(i).NE.0)) THEN
          WRITE(13,111) -nx(i),-ny(i),nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,111) -nx(i),ny(i),-nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,111) nx(i),-ny(i),-nz(i)
       ENDIF
    ENDDO

    DO i=1,npt
       IF((nx(i).NE.0).AND.(ny(i).NE.0).AND.(nz(i).NE.0)) THEN
          WRITE(13,111) -nx(i),-ny(i),-nz(i)
       ENDIF
    ENDDO

111 FORMAT(5(E20.10,1x))
112 FORMAT(5(E20.10,1x))

    CLOSE(11, status='delete')
    CLOSE(12, status='delete')
    CLOSE(14, status='delete')
    CLOSE(16, status='keep')

    !DEALLOCATE(nx, ny, nz, c1)

  END SUBROUTINE generrpr



  SUBROUTINE dcgen(nx,ny,nz,cd)

    IMPLICIT NONE

    REAL(prcn), INTENT(in) ::  nx,ny,nz
    REAL(prcn), INTENT(out):: cd(3)
    REAL(prcn) :: r 

    r=dsqrt(DBLE(nx*nx+ny*ny+nz*nz))

    cd(1)=DBLE(nx)/r
    cd(2)=DBLE(ny)/r
    cd(3)=DBLE(nz)/r
  END SUBROUTINE dcgen



END MODULE geom_init
  
