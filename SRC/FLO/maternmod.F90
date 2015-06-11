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


MODULE maternmod
#include "ibm.h"
  USE precision 
  USE constants 
  Use collision_mod
  USE randomno, ONLY : uni_dist
!	use mypost_process
	use postproc_funcs
  IMPLICIT NONE 
  PRIVATE 
  ! REAL(prcn) :: epan
  PUBLIC ::  matern,matern_check_nd
CONTAINS 
  SUBROUTINE matern(ndim,iborddel,L,numdens_p,r_0,hc_fac,nsim,npmax,xc,np) !,final_vol_frac)


!!$c     physical parameters
!!$
!!$c     numdens_P: expected number density of Poisson process
!!$c     numdens_M: expected number density of Matern process
!!$c     numdens_est: estimated number density of Matern process
!!$c     numdens_sq_est: estimated second moment of Matern process

    IMPLICIT NONE
    INTEGER, INTENT(in) :: ndim,npmax,nsim
    LOGICAL, INTENT(IN) :: iborddel

    INTEGER, PARAMETER :: ndim0=3, nrmax = 350, nsimmax=1000!ndim0
    ! set to two ... change for 3d case 
    REAL(prcn),INTENT(in) :: numdens_P, L(ndim),r_0, hc_fac
    REAL(prcn),INTENT(out) , DIMENSION(:,:)::  xc!(npmax,ndim0)
!	real(prcn), intent(out) :: final_vol_frac
    INTEGER, INTENT(out) :: np 
    REAL(prcn) ::  xp(npmax,ndim0)
    REAL(prcn) :: numdens_M, palm_p, numdens_est_isim, numdens_est,&
         numdens_sq_est_isim, numdens_sq_est
!!$
    REAL(prcn), DIMENSION(npmax) :: wt, rp, frand,&
         & e, mark, wttemp,& 
         & frandr, dist  
    REAL(prcn) :: four_thirdpi
    INTEGER  :: i, idim, ir, j, irmin, irmax,  npmatern,&
         isim, npairs_tooclose, ndel, nnewdel, ndelbord,&
         iunit,k,icount,icnt,ibin, &
         incr,kpoint,iflag,npairs,iscount, nptot, NPARTICLES 
    INTEGER, DIMENSION(npmax)::ifpon
    Real(prcn), DIMENSION(npmax,ndim0) ::sp
    REAL(prcn) ::  sumy,&
         &c, h, rho_est_isim(nrmax), dij,&
         & Adenom, W1min, Vol, & 
         Lint(ndim0), r, temp, rmin, rmax, dr, &
         rmin_contrib, rmax_contrib, rball, vol0, rvol0, &
         rho_est(nrmax), Vol_trunc, H0(ndim0,ndim0)

    REAL(prcn) ::  sum_coord,rho_est_md(nrmax), &
         rho_est_md_isim(nrmax),vol_ir,dij2, &
         mindist,xp_temp(ndim0),r2,r1, &
         gamma(nrmax),Ugamma,rho_analy(nrmax), gofr_avg(nrmax),&
         & rho_est_avg(nrmax), dxeff
    integer :: ip

	real(prcn), allocatable :: radtmp(:)
	real(prcn) :: tmp_vec(ndim), dist1, min_dist
	
	logical, allocatable :: contact(:,:)
	real(prcn) :: max_overlap

    !define constants

    !open(10,file='matin.data',status='old')

    OPEN(2,file='xp_poisson.dat',form='formatted')
    OPEN(3,file='rhoest.dat',form='formatted')
    OPEN(14,file='xp_matern.dat',form='formatted')
    OPEN(12,file='xp_matern.cfg',form='formatted')
    OPEN(20,file='cluster.dat',form='formatted')

    cgrid_fac = 1.5d0
    MN = 20
    FACTOR_RLM = HC_FAC
    PARTICLES_FACTOR = 1.0
    
    dxeff = two*r_0*cgrid_fac

    cy = MAX(NINT(L(2)/dxeff),2)
    cx = MAX(NINT(L(1)/dxeff),2)
    cz = MAX(NINT(L(3)/dxeff),2)

!!$    cx = NINT(L(1)/(cgrid_fac*two*r_0))
!!$    cy = NINT(L(2)/(cgrid_fac*two*r_0))
!!$    cz = NINT(L(3)/(cgrid_fac*two*r_0))

    ALLOCATE(cgrid(cx,cy,cz,3))
    dxc = L(1)/REAL(cx-1,prcn)
    dyc = L(2)/REAL(cy-1,prcn)
    dzc = L(3)/REAL(cz-1,prcn)

    DO k = 1, cz
       do j = 1, cy
          do i = 1, cx 
             cgrid(i,j,k,1) = (i-1)*dxc
             cgrid(i,j,k,2) = (j-1)*dyc
             cgrid(i,j,k,3) = (k-1)*dyc
          end do
       end do
    end DO
!!$    PRINT*,'doml ', L(1), L(2)
!!$    PRINT*,'cx, cy, cz = ', cx,cy,cz
!!$    PRINT*,'dxc, dyc, dzc = ', dxc,dyc,dzc



    IMAX = CX - 1
    JMAX = CY - 1
    KMAX = CZ - 1

    IMAX1 = IMAX + 1
    JMAX1 = JMAX + 1
    KMAX1 = KMAX + 1

    IMAX2 = IMAX1+1
    JMAX2 = JMAX1+1
    KMAX2 = KMAX1+1

    IMAX3 = IMAX2
    JMAX3 = JMAX2
    KMAX3 = KMAX2

    IMIN1 = 2
    JMIN1 = 2
    KMIN1 = 2

    IMIN2 = 1
    JMIN2 = 1
    KMIN2 = 1

!!$    PRINT*,'IMAX, JMAX, KMAX = ', IMAX, JMAX, KMAX
!!$    PRINT*,'IMIN1, JMIN1, KMIN1 = ', IMIN1, JMIN1, KMIN1
!!$    PRINT*,'IMIN2, JMIN2, KMIN2 = ', IMIN2, JMIN2, KMIN2
!!$
!!$    PRINT*,'IMAX1, JMAX1, KMAX1 = ', IMAX1, JMAX1, KMAX1
!!$    PRINT*,'IMAX2, JMAX2, KMAX2 = ', IMAX2, JMAX2, KMAX2
!!$    PRINT*,'IMAX3, JMAX3, KMAX3 = ', IMAX3, JMAX3, KMAX3

    ALLOCATE(XE(IMAX2), YN(JMAX2), ZT(KMAX2))

    XE = ZERO
    YN = ZERO
    ZT = ZERO
    XE(1) = ZERO
    YN(1) = ZERO
    ZT(1) = ZERO

    DO I = IMIN1, IMAX2
       XE(I) = XE(I-1) + DXC
       !PRINT*,'XE = ', I, XE(I)
    END DO

    DO J  = JMIN1, JMAX2
       YN(J) = YN(J-1) + DYC
    END DO

    DO K = KMIN1, KMAX2
       ZT(K) = ZT(K-1) + DZC
    END DO
    MAXNEIGHBORS = MN + 6

    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A40,2(2x,g17.8))') 'XLENGTH CALC =', XE(IMAX1) - XE(IMIN2)
       WRITE(*,'(A40,2(2x,g17.8))') 'YLENGTH CALC =', YN(JMAX1) - YN(JMIN2)
       WRITE(*,'(A40,2(2x,g17.8))') 'ZLENGTH CALC =', ZT(KMAX1) - ZT(KMIN2)
    end if

    !OPEN(1,file='seed.d')
    four_thirdpi=4.0*pi/3.0
    c = 0.15
    rmax = 0.0

    
    Vol=1.0
    DO idim = 1, ndim
       Vol = Vol*L(idim)
       rmax = rmax + L(idim)**2
    ENDDO
    rmax = SQRT(rmax)         !maximum distance inside the domain
    rmin = 0.0                !lower limit of g(r) graph
    dr = (rmax-rmin)/float(nrmax) !step in r
    DO ir = 1, nrmax
       rho_est(ir) = 0.0      !initialize
       rho_est_md(ir) = 0.0      !initialize
    ENDDO


    rball = hc_fac*r_0            !hard core distance

    IF(ndim.EQ.3)THEN
       vol0 = four_thirdpi*(rball**3) 
    ELSEIF (ndim.EQ.2)THEN
       vol0 = pi*rball**2     
    ELSEIF (ndim.EQ.1)THEN 
       vol0 = 2.*rball 
    ENDIF
    ! write(*,*)'Check: ndim = ',ndim

    !Calculate number density of matern process and palm probability
    numdens_M = (one-EXP((-1.)*numdens_P*vol0))/vol0
    palm_p    = (one-EXP((-1.)*numdens_P*vol0))/(vol0*numdens_P)

	write (*,*) "Expected number density of Matern process = ", numdens_M
	write (*,*) "Expected number density of Palm process   = ", palm_p


    !/ distribute particles according to specified point process
    !choose Poisson process to begin with

    !define domain to be unit square or unit cube

    !initialize particle properties

    !According to SKM p. 29 sequence of random points uniformly
    !distributed in the hypercube [0,1]^d is produced by
    !x_i = (z_(i-1)*d+1, ..., z_i*d ) for i=1,2,...

    !For a plane this means: generate z_n a sequence uniformly
    !distributed in [0,1] (ranu2 will do this)
    !x_i = (z_2*i-1, z_2*i) for i=1,2,...
    !i.e., x_1 = (z1, z2), x_2 = (z3,z4) etc
    !For 3D
    !x_1 = (z1, z2, z3), x_2 = (z4, z5, z6) etc

    iscount = 1

    numdens_est=0.0           !estimated number density
    numdens_sq_est=0.0        !estimated number density squared

!	final_vol_frac = zero

    DO isim=1, nsim           !begin loop for simulations

       DO ir = 1, nrmax
          rho_est_isim(ir) = 0.0 
          rho_est_md_isim(ir) = 0.0
       ENDDO

       CALL uni_dist(frand)

       !call ranu2(frand,1,npmax,1,1,npmax)

       DO i = 1, npmax
          ifpon(i) = 1
       ENDDO

       sumy = 0.0
       i=0
       DO WHILE (sumy .LE. numdens_P*Vol-1) 
          i = i + 1
          IF(i.GT.npmax) then
             PRINT*,'Error in matern....Need to increase the value of npmax.... Current npmax = ', npmax,' i = ',i 
             Write(*,*)'SUMY = ', SUMY, 'NUM : ', numdens_P*vol
          END IF
          e(i) = -LOG(frand(i)) !generate exponential random variate
          sumy = sumy + e(i)  !add each random variate
          ! print*,'i=',i,'e(i) =',e(i),'sumy =',sumy, &
          !      'numdens=',numdens_P, 'ndvol=',numdens_p*vol
       ENDDO
       Write(*,*)'SUMY = :', sumy, 'Num : ', numdens_P*vol, 'NPMAX: ', npmax
       np = i                 !Poisson distributed random variate
       write(*,*)'Value of np..',np,'..Simulation...',isim

       IF ( np .GT. npmax) THEN
          WRITE(*,*)'ERROR: np > npmax'
          WRITE(*,*)'       Change parameter npmax to >=',np
          STOP
       ENDIF

       DO idim=1,ndim
          CALL uni_dist(xp(1:np,idim))
       ENDDO !generate np uniformly 
       !random variates  
       IF(ndim.LT.ndim0) xp(:,ndim+1:ndim0) = zero !RG 09/04/06 

       !wt - particle statistical weight (number weighting here)

       IF (isim .EQ. 1) THEN
          !WRITE(2,*)'@ title "Particle positions and radius for Pois&
           !    &son distribution"' 
          ! WRITE(3,*)'@ title "Estimate of second order stats for Mat&
          !      &ern process:"' 
          ! WRITE(3,*) '# r, rho(r), g(r)'
          !WRITE(4,*)'@ title"Particle positions and radius for Mater&
          !     &n distribution"' 
       ENDIF

       DO i = 1, np
          rp(i) = r_0
          wt(i) = 1.0/float(np) !equal weights to all
          DO idim=1, ndim
             xp(i,idim) = L(idim)*xp(i,idim) !positions of particles
             !uniformly distributed
          ENDDO
       ENDDO

       !Added later from matern_new.f RG 09/04/06
       !GOTO 100

       h = c/SQRT(numdens_M)
!!$       Print*,' h , c, Vol = ', h, c, Vol
!!$       Print*, 'lam_m, lam_p = ',numdens_m, numdens_p
!!$       Print*, 'r_o, rball, vol0 = ', r_0, rball, vol0, isim
       !Read(*,*)
       !loop over all particle pairs (excluding i=j)


       CALL GRID_BASED_NEIGHBOR_SEARCH_MAT(NP,XP(1:NP,1:3), r_0, IFPON(1:NP), FRAND(1:NP))

       npairs_tooclose=0
       ndel=0
       nnewdel=0
       
       DO i = 1, np
          IF(IFPON(I).EQ.0) THEN 
             nnewdel=nnewdel+1
             ndel=ndel+1
          end IF
       end DO
!!$       
!!$
!!$       DO i = 1, np
!!$          DO j = 1, np
!!$             IF ( j.NE.i) THEN
!!$                ! May 4, 2004, SS: See note on p 165 SKM: "Here
!!$                ! care is necessary: 
!!$                !points thinned out can nevertheless thin out other
!!$                ! points with greater marks". Therefore ifpon
!!$                ! should not be used in this loop, and we need to
!!$                ! loop over all distinct pairs. calculate the inter
!!$                !-point distance dij 
!!$                dij = 0.0
!!$                DO idim = 1, ndim
!!$                   dij = dij + (xp(i,idim) - xp(j,idim))**2
!!$                ENDDO
!!$                dij = SQRT(dij)
!!$                IF (dij .LE. rball) THEN
!!$                   npairs_tooclose=npairs_tooclose+1
!!$                   IF ( frand(j) .LT. frand(i) ) THEN
!!$                      IF (ifpon(i) .EQ. 1) THEN
!!$                         nnewdel=nnewdel+1
!!$                      ENDIF
!!$                      ifpon(i) = 0
!!$                      ndel=ndel+1
!!$                   ENDIF
!!$                ENDIF
!!$             ENDIF
!!$          ENDDO
!!$          !Print*, 'dij = ', dij, rball
!!$          
!!$       ENDDO
       ndelbord=0

       IF(iborddel) THEN 
          PRINT*,'removing all points in edge border 2r_0 of cube'

          DO i = 1, np
             DO idim=1, ndim !Change RG 09/04/06
                IF ( (ABS(xp(i,idim)) .LE. rp(i)) .OR. (ABS(xp(i,idim)-L(idim)) .LE. rp(i)) ) THEN 
                   ifpon(i) = 0
                   ndelbord=ndelbord+1
                ENDIF
             ENDDO
          ENDDO
       end IF

       !separate calculation of rho_2
!!$       DO i = 1, np
!!$          DO j = 1, np
!!$             IF ( j.NE.i .AND.(ifpon(i).EQ.1) .AND. (ifpon(j).EQ.1)&
!!$                  & ) THEN 
!!$                !calculate the inter-point distance dij
!!$                dij = 0.0
!!$                Adenom = 1.0
!!$                DO idim = 1, ndim
!!$                   dij = dij + (xp(i,idim) - xp(j,idim))**2
!!$
                   !calculate contributions to the denominator
                   ! A(Wxj \int Wxi)  our window is a unit cube 

                   !calculate the intersection of the segments
                   ! (xp(i,idim), xp(i,idim) + L(idim)) and (xp(j
                   !,idim), xp(j,idim)+L(idim)) 

                   !Find global min of Wi and Wj (in x,y,z succ),
                   ! in x call this xmin; if xp(i, ) = xmin we're
                   ! in i, else in j; call this window 1 find the
                   ! diff betn max window 1 - min window 2 if this
                   ! diff is positive, there is a nonzero
                   ! intersection, else this diff is negative and
                   ! the windows don't intersect 

!!$                   W1min = MIN(xp(i,idim),xp(j,idim))
!!$                   IF ( xp(i,idim) .EQ. W1min ) THEN
!!$                      !the window 1 is the i window
!!$                      !Lint is the length of the intersecting
!!$                      ! segment in dim idim 
!!$                      !L(idim) = 1.0 in our case
!!$                      Lint(idim) = MAX(zero,(xp(i,idim)+L(idim)-xp(j&
!!$                           &,idim))) 
!!$                   ELSE
!!$                      !the window 1 is the j window
!!$                      Lint(idim) = MAX(zero,(xp(j,idim)+L(idim)-xp(i,idim)))
!!$                   ENDIF
!!$                   Adenom = Adenom*Lint(idim)
!!$                ENDDO
!!$                dij = SQRT(dij)
!!$                Print*,'xp = ', xp(i,1), xp(j,1)
!!$                Print*,'xp = ', xp(i,2), xp(j,2)
!!$                Print*, 'Adenom = ', i,j, Adenom, lint(1), lint(2)

                !Use the Epanecnikov kernel; note that each
                ! interpoint distance contributes to r bins that
                ! are within the kernel width +/- h 

                !kernel width parameter h is recommended as c/sqrt(lambda)
                !take c = 0.15

                !calculate rmin and rmax that this pair contributes to

!!$                rmin_contrib = dij - h
!!$                rmax_contrib = dij + h

                !for each bin r(ir) such that rmin <= r(ir) <=rmax
                ! calculate the value of rho_est 

!!$                irmin = int((rmin_contrib-rmin)/dr)
!!$                irmax = int((rmax_contrib-rmin)/dr )
!!$                irmin = (rmin_contrib-rmin)/dr 
!!$                irmax = (rmax_contrib-rmin)/dr 
!!$                    write(*,*)'Pair i=',i,'j = ',j, dij,h
!!$                    write(*,*)'Adenom = ',Adenom
!!$                    write(*,*)'rmin_contrib = ',rmin_contrib
!!$                    write(*,*)'rmax_contrib = ',rmax_contrib
                !Print*,'irmin, irmax  = ', irmin, irmax
!!$                DO ir = irmin, irmax
!!$                   r = rmin + ir*dr
!!$                   temp = r - dij
!!$                   rho_est_isim(ir) = rho_est_isim(ir) + epan(temp&
!!$                        &,h)/Adenom 
!!$                     write(*,*)'temp = ',temp,'h = ',h,' contrib='&
!!$                          &,epan(temp,h)/(Adenom)
!!$                     write(*,*)'ir = ',ir,' rho_est(ir) = ',rho_est(ir)
!!$                ENDDO
!!$             ENDIF
!!$          ENDDO
!!$       ENDDO

!!$       DO ir = 1, nrmax
!!$          r = rmin + ir*dr
!!$          if(ndim.eq.1)then
!!$             rho_est_isim(ir) = rho_est_isim(ir)/(2.*r)
!!$          elseif(ndim.eq.2)then
!!$             rho_est_isim(ir) = rho_est_isim(ir)/(2.*pi*r) 
!!$          elseif(ndim.eq.3)then
!!$             rho_est_isim(ir) = rho_est_isim(ir)/(4.*pi*r**2)
!!$          endif
!!$          !rho_est_isim(ir) = rho_est_isim(ir)/(4.*pi*r**2)
!!$          rho_est(ir) = rho_est(ir) + rho_est_isim(ir)/float(nsim)
!!$          WRITE(30,'(10(2x,e20.10))')r, rho_est(ir), rho_est(ir)/(numdens_M**2),rho_est_isim(ir)  
!!$       ENDDO

       npmatern = 0
       DO i = 1, np
          IF (ifpon(i) .EQ. 1) THEN
             npmatern = npmatern + 1
             xc(npmatern,1:ndim) = xp(i,1:ndim)
             !WRITE(4,*)(xp(i,idim),idim=1,ndim),rp(i)
          ENDIF
       ENDDO

       !total number of particles ! entered by Madhu - for plotting
       ! purposes 

       nptot = nptot + npmatern

       iunit = iunit + 1

       !2/ compute particle number density mean

       Vol_trunc=1.0
       DO idim = 1, ndim
          Vol_trunc = Vol_trunc*(L(idim)-0.0*rball)
       ENDDO
       numdens_est_isim = float(npmatern)/Vol_trunc
       numdens_sq_est_isim =float(npmatern)*float(npmatern-1)/Vol_trunc**2 
       numdens_est = numdens_est + numdens_est_isim/float(nsim)
       numdens_sq_est = numdens_sq_est +numdens_sq_est_isim/float(nsim) 
!!$
!!$       IF (isim .EQ. 1) THEN
!!$          rvol0 = 1./vol0
!!$          WRITE(6,101)numdens_P, numdens_M, palm_p, r_0, rball&
!!$               &,(L(idim),idim=1,ndim), vol0, rvol0 
!!$          PRINT*, 'ndim  = ', ndim
!!$101       FORMAT(1h ,/2x,' Particle Distribution by Point Process Si&
!!$               &mulation',/,2x,' Matern Process',/,5x,' Physical Par&
!!$               &ameters',/, 8x,' Poisson Number density',t60,f8.3,/&
!!$               & 8x,' Matern Number density',t60,f8.3,/8x,' Palm ret&
!!$               &aining probability',t60,f8.3,/8x,' Radius of particl&
!!$               &es',t60,f8.3,/8x,' Hard core distance: h=',t60,f8.3&
!!$               &,/8x,' Box dimensions',t44,3f8.3,/8x,' Constant c=b_&
!!$               &d h^d',t60,e10.3,/8x,' Max n_m for this h=1/c=1/b_d &
!!$               &h^d',t60,e10.3) 
!!$       ENDIF
!!$       WRITE(6,102)isim, numdens_est_isim, np, npmatern&
!!$            &,float(npmatern)/float(np), npairs_tooclose/2, ndel,&
!!$            & nnewdel, ndelbord 
!!$102    FORMAT(1h ,4x,' Simulation number',t60,i8,/, 8x,' Estimated Ma&
!!$            &tern Number density',t60,f8.3,/ 8x,' Number of Poisson &
!!$            &particles:',t60,i8,/8x,' Number of Matern particles:'&
!!$            &,t60,i8,/8x,' Estimate of Palm retaining probability'&
!!$            &,t60,f8.5,/8x,' Pairs too close:',t60,i8,/8x,' Total nu&
!!$            &mber of deletions:',t60,i8,/8x,' Number of new deletion&
!!$            &s:',t60,i8,/8x,' Number of border deletions:',t60,i8) 
!!$
!!$       !3/ plot
!!$
!!$    ENDDO !isim 

!!$    WRITE(13,*)nptot
!!$
!!$    WRITE(6,103)numdens_est
!!$103 FORMAT(1h ,/ 8x,' MIS Estimated Matern Number density',t60,f8.3&
!!$         &,/) 

!!$    rho_est_avg = zero 
!!$    gofr_avg = zero 
!!$    
!!$    do j=1,nrbins
!!$       rho_est_avg(j) = rho_est_avg(j) + 1./dble(nsim)*sum(rho_est(1:nsim,j))
!!$       gofr_avg(j) = gofr_avg(j) + 1./dble(nsim)*sum(gofr(1:nsim,j))
!!$    end do


!!$    DO ir = 1, nrmax
!!$       r = rmin + ir*dr
!!$       WRITE(3,'(10(2x,e20.10))')r/rball, rho_est(ir), rho_est(ir)&
!!$            &/(numdens_sq_est),rho_est_isim(ir)/numdens_sq_est&
!!$            &,rho_est_isim(ir)  



!write (*,*) xc
!read (*,*)

		CALL screen_separator(30,'^')
		write (*,*) "NUMBER OF PARTICLES = ", npmatern
          write (*,*) "MATERN VOLUME FRACTION = ", npmatern*fourthirdpi*(r_0**3)/(l(1)*l(2)*l(3))

!		final_vol_frac = final_vol_frac + npmatern*fourthirdpi*(r_0**3)/(l(1)*l(2)*l(3))

		if (.not.allocated(contact)) allocate(contact(nbody,nbody))
		!call calculate_gofr_homog(npmatern,xc(1:npmatern,1:ndim), contact, my, mxf, nbins, .true., gofr_mis(isim,1:nbins), rad_bin(1:nbins), max_overlap)


		if (.not.allocated(radtmp)) allocate(radtmp(npmatern))
		radtmp = r_0

!		call interstitial_dist(xc(1:npmatern,1:ndim), radtmp(1:npmatern), int_dist(isim), npmatern)
		if (allocated(radtmp)) deallocate(radtmp)


		min_dist = 1000.
		do i=1, npmatern-1
			do j=i+1, npmatern
				tmp_vec(:) = xc(i,:)-xc(j,:)
				dist1 = sqrt(dot_product(tmp_vec,tmp_vec))

				if (dist1<min_dist) min_dist = dist1
			enddo
		enddo
		write (*,*) "MINIMUM SEPARATION = ", min_dist


		CALL screen_separator(30,'-')
    ENDDO

    np = npmatern
!	final_vol_frac = final_vol_frac/nsim

    DEALLOCATE(CGRID, XE, YN, ZT)
100 CONTINUE
  END SUBROUTINE matern


  subroutine matern_check_nd

    USE precision 
    USE constants
    implicit none 
    !This subroutiune is just a test routine to compare the Metern
    ! distribution for the case of 1-d, 2-d, and 3 -d, distributions. 
    !It was observed  that for 1-d case, the metern density was very
    ! less compared to input Poisson number density for same values of
    ! hard core distance, particle radius, and poisson number density
    !RG 09/21/06

    Real(prcn) :: hc, rp, vol0, lambdap, lambdam 

    Real(prcn) :: h1d, h2d, h3d, vol1d, vol2d, vol3d 
    Real(prcn) :: rhat1d, rhat2d, rhat3d, rhat, hhat, c(3), thetahat(3)&
         &, dhhat, drhat, lambmbylambp(3), rhatnew

    Integer :: i,j,k , nhhat, nrhat 
    nhhat = 20
    nrhat = 20
    dhhat = one/(nhhat-1)
    drhat = (one-0.01)/(nrhat-1)
    OPEN(201,file='volfr.dat',form='formatted')
    write(201,*)'VARIABLES= ',' "rhat" ',' "hhat" ',' "theta1" ',' "thet&
         &a2" ',' "theta3" ' 
    write(201,*) 'ZONE F=POINT, I=',nrhat,  ', J=',nhhat

    OPEN(202,file='lambda_ratio.dat',form='formatted')
    write(202,*)'VARIABLES= ',' "rhat" ',' "hhat" ',' "nmbynp1" ',' "nmb&
         &ynp2" ',' "nmbynp3" '  
    write(202,*) 'ZONE F=POINT, I=',nrhat,  ', J=',nhhat
!!$    do k = 1,nk

    c(3) = (four*pi)/three
    c(2) =  pi
    c(1) = two
    rp = 0.007d0
    do hhat = 2,3 , dhhat
       do rhat = 0.01,1, drhat

          do i = 1,3 
             rhatnew = ((rhat/rp)**(3/i))*rp
             thetahat(i) =  (one-EXP((-1.)*c(i)*(hhat**i)*(rhatnew**i)))/hhat&
                  &**i
             lambmbylambp(i)=(one-EXP((-1.)*c(i)*(hhat**i)*(rhatnew**i)))&
                  &/(c(i)*(hhat**i)*(rhatnew**i))
          end do

          write(201,'(10(1x,E12.5))') rhat, hhat, thetahat(1),&
               & thetahat(2), thetahat(3)
          write(202,'(10(1x,E12.5))') rhat, hhat,lambmbylambp(1)&
               &,lambmbylambp(2),lambmbylambp(3)
       end do

    end do

    close(201)
    close(202)
  end subroutine matern_check_nd

  SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH_MAT(PARTICLES,DES_POS_NEW, RAD1,IFPON, FRAND)

    IMPLICIT NONE
    LOGICAL PER_COND, ALREADY_NEIGHBOURS
    INTEGER I, II, LL, CO, NI, TEMP, JJ, KK , J, K, NEIGH_L, PARTICLES, DIMN, L_MAX, PNO_MAX
    INTEGER KM1, KP1, IM1, IP1, JM1, JP1, PNO, NPG, PC(3), IP2, NLIM
    REAL(PRCN), INTENT(IN) :: DES_POS_NEW(PARTICLES,3), RAD1, FRAND(PARTICLES)
    DOUBLE PRECISION  DIST(NDIM), DISTMAG, R_LM, LX, LY, LZ, XPER_FAC, YPER_FAC, ZPER_FAC,  CORD_PNO(NDIM), MAX_OVERLAP, TMP_OVERLAP,  DES_RADIUS(PARTICLES)
    INTEGER , INTENT(OUT) :: IFPON(PARTICLES)
    real(prcn) :: tempx, tempy, tempz
    DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

    INTEGER ::  npic ,pos, IP

    INTEGER, DIMENSION(1:cx,1:cy,1:cz):: iicount
    
!!$    PRINT*, 'IN CELL LINKED LIST SEARCH', INTX_PER, INTY_PER, INTZ_PER, PARTICLES
!!$    PRINT*,'IMAX, JMAX, KMAX = ', IMAX, JMAX, KMAX
!!$    PRINT*,'IMIN1, JMIN1, KMIN1 = ', IMIN1, JMIN1, KMIN1
!!$    PRINT*,'IMIN2, JMIN2, KMIN2 = ', IMIN2, JMIN2, KMIN2
!!$
!!$    PRINT*,'IMAX1, JMAX1, KMAX1 = ', IMAX1, JMAX1, KMAX1
!!$    PRINT*,'IMAX2, JMAX2, KMAX2 = ', IMAX2, JMAX2, KMAX2
!!$    PRINT*,'IMAX3, JMAX3, KMAX3 = ', IMAX3, JMAX3, KMAX3
!!$
    
    
    ALLOCATE(pic(IMAX2,JMAX2,KMAX2))
    ALLOCATE(cnd(IMAX2,JMAX2,KMAX2))

    ALLOCATE(PIJK(PARTICLES,3))

    Allocate(  NEIGHBOURS (PARTICLES, MAXNEIGHBORS) )

    
    
    DO  k  = 1,KMAX2!MAX(KMAX1-1,1)
       DO  j  = 1,JMAX2
          DO  i  = 1,IMAX2
             NULLIFY(pic(i,j,k)%p)

             cnd(i,j,k) = 0
          end DO
       end DO
    end DO


    DO IP = 1, PARTICLES
       tempx = (des_pos_new(ip,1))
       tempy = (des_pos_new(ip,2))
       tempz = (des_pos_new(ip,3))
       PIJK(IP,1) = INT(tempx/dxc)+2
       PIJK(IP,2) = INT(tempy/dyc)+2
       PIJK(IP,3) = INT(tempz/dzc)+2
       !PRINT*,'tempx = ', PIJK(IP,:)

       cnd(PIJK(IP,1),PIJK(IP,2),PIJK(IP,3)) = cnd(PIJK(IP,1),PIJK(IP,2),PIJK(IP,3)) + 1

    end DO


    DO  k = 1,KMAX2        !MAX(KMAX1-1,1)
       DO  j = 1,JMAX2
          DO  i = 1,IMAX2
             NPIC = CND(i,j,k)
             !PRINT*,'NPIC = ', NPIC
             IF (ASSOCIATED(pic(i,j,k)%p)) THEN
                IF (npic.NE.SIZE(pic(i,j,k)%p)) THEN
                   DEALLOCATE(pic(i,j,k)%p)
                   IF (npic.GT.0) ALLOCATE(pic(i,j,k)%p(npic))
                ENDIF
             ELSE
                IF (npic.GT.0) ALLOCATE(pic(i,j,k)%p(npic))
             ENDIF

          end DO
       end DO
    end DO

    iicount(:,:,:) = 1

    DO ip = 1, PARTICLES
       PC(:) =  PIJK(IP,1:3)
       pos = iicount(pc(1),pc(2),pc(3))
       pic(pc(1),pc(2),pc(3))%p(pos) = ip
       !PRINT*,'pic = ', pic(pc(1),pc(2),pc(3))%p(pos)
       iicount(pc(1),pc(2),pc(3)) = &
            & iicount(pc(1),pc(2),pc(3)) + 1
    ENDDO

    DO I = 1, PARTICLES
       IFPON(I) = 1
       DES_RADIUS(I) = RAD1
       !PRINT*,'rad = ', DES_RADIUS(I)
       DO II = 1, MAXNEIGHBORS
          NEIGHBOURS(I,II) = -1
       END DO
       NEIGHBOURS(I,1) = 0
    END DO

    MAX_RADIUS = RAD1
    LX = XE(IMAX1) - XE(1)
    LY = YN(JMAX1) - YN(1)
    LZ = ZT(KMAX1) - ZT(1)
    !PARTICLES = NBODY
    DIMN = NDIM
    MAX_OVERLAP = SMALL_NUMBER
    PNO_MAX = 0
    L_MAX= 0

    !WRITE(*,'(A,3(2x,g12.5,/))') ' DOMAIN LENGTH = ', LX, LY, LZ
    DO LL = 1, PARTICLES
       IF(IFPON(LL).EQ.0) GOTO 999

       PC(:) =  PIJK(LL,1:3)
       !PRINT*,'PC = ', PC
       II = PC(1)
       JJ = PC(2)
       KK = PC(3)

       IP1 =II+1! MIN(IMAX1, II+1)
       IM1 =II-1! MAX(IMIN1, II-1)
       JP1 =JJ+1! MIN(JMAX1, JJ+1)
       JM1 =JJ-1! MAX(JMIN1, JJ-1)
         
       IF(DIMN.EQ.3) THEN 
          KP1 =KK+1! MIN(KMAX1, KK+1)
          KM1 =KK-1! MAX(KMIN1, KK-1)
       end IF

       DO KK = KM1, KP1
          DO JJ = JM1, JP1
             DO II = IM1, IP1

                I = II
                J = JJ
                K = KK

                XPER_FAC = 0
                YPER_FAC = 0
                ZPER_FAC = 0
                PER_COND = .FALSE.
                IF(II.GT.IMAX1) THEN 
                   IF(INTX_PER) THEN 
                      I = IMIN1
                      XPER_FAC = one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true EAST',I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      I = IMAX1
                   ENDIF
                ENDIF

                IF(II.LT.IMIN1) THEN 
                   IF(INTX_PER) THEN 
                      I = IMAX1
                      XPER_FAC = -one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true WEST', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE 
                      I = IMIN1
                   ENDIF
                ENDIF

                IF(JJ.GT.JMAX1) THEN 
                   IF(INTY_PER) THEN 
                      J = JMIN1
                      YPER_FAC = one
                      PER_COND = .true.
                     ! WRITE(*,*) 'cond true NORTH', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      J = JMAX1
                   ENDIF
                ENDIF

                IF(JJ.LT.JMIN1) THEN 
                   IF(INTY_PER) THEN 
                      J = JMAX1
                      YPER_FAC = -one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true SOUTH', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      J = JMIN1
                   ENDIF
                ENDIF

                IF(DIMN.EQ.3) THEN 
                   IF(KK.GT.KMAX1) THEN 
                      IF(INTZ_PER) THEN
                         K = KMIN1
                         ZPER_FAC = one
                      ELSE
                         K = KMAX1
                      ENDIF
                   ENDIF

                   IF(KK.LT.KMIN1) THEN 
                      IF(INTZ_PER) THEN
                         K = KMAX1
                         ZPER_FAC = -one
                      ELSE
                         K = KMIN1
                      ENDIF
                   ENDIF
                ENDIF

                If (ASSOCIATED(PIC(I,J,K)%p)) then
                   NPG = SIZE(PIC(I,J,K)%p)
                   !PRINT*, 'LL = ', LL, I,J,K, NPG
                   !DO IP2 = 1, NPG
                   !   PRINT*, 'PNO2 = ', PIC(I,J,K)%p(ip2)
                   !end DO

                Else
                   NPG = 0
                Endif
                
                Do IP2 = 1,NPG
                   PNO = PIC(I,J,K)%p(ip2)
                   

                   !  PRINT*,'LL, PNO = ', LL, PNO, NPG, I,J,K
                   R_LM = DES_RADIUS(LL)! + DES_RADIUS(PNO)
                   R_LM = FACTOR_RLM*R_LM
                   CORD_PNO(1) = DES_POS_NEW(PNO,1) + XPER_FAC*(LX)
                   CORD_PNO(2) = DES_POS_NEW(PNO,2) + YPER_FAC*(LY)!+ONE)
                   IF(DIMN.EQ.3) THEN 
                      CORD_PNO(3) = DES_POS_NEW(PNO,3) + ZPER_FAC*(LZ)!+ONE)
                   ENDIF

                   
                   DIST(:) = CORD_PNO(:) - DES_POS_NEW(LL,:)
                   DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
                   
                   ALREADY_NEIGHBOURS = .FALSE.
                   
                   !IF(LL.EQ.1.AND.PNO.EQ.2) THEN 
                   !   WRITE(*,*)'CORD=', CORD_PNO(1)
                   !   WRITE(*,*)'DISTMAG', DISTMAG, R_LM-DISTMAG
                   !ENDIF
                   

                   IF(R_LM - DISTMAG.gt.SMALL_NUMBER.AND.(.NOT.ALREADY_NEIGHBOURS)) THEN 
                      
                      TMP_OVERLAP = ((DES_RADIUS(LL) + DES_RADIUS(PNO))-DISTMAG)/(DES_RADIUS(LL) + DES_RADIUS(PNO))
                      TMP_OVERLAP = TMP_OVERLAP*100
                      IF(TMP_OVERLAP.GT.MAX_OVERLAP) THEN 
                         
                         MAX_OVERLAP = MAX(MAX_OVERLAP, TMP_OVERLAP)
                         L_MAX = LL
                         PNO_MAX = PNO
                         
                      end IF
                      


!!$                           IF(PER_COND) THEN 
!!$                              WRITE(*,*) 'pC = ', pc
!!$                              WRITE(*,*) 'II, JJ = ', II, JJ
!!$                              WRITE(*,*) 'I, J = ', I, J
!!$                              WRITE(*,*) 'XYPER_FAC ', XPER_FAC, YPER_FAC
!!$                              WRITE(*,*) 'DES_VEL_NEW = ', DES_POS_NEW(PNO,:)
!!$                              WRITE(*,*) 'MODIFIED POSITION = ', CORD_PNO(:)
!!$                           ENDIF



                         
                      IF(FRAND(PNO).LT.FRAND(LL)) THEN 
                         
                         ifpon(LL) = 0
!!$                         IF(PER_COND)THEN 
!!$                            PRINT*,'PART DEACTIVATED', PER_COND
!!$                         end IF
                         
                         goto 999
                      end IF
                      
                   end IF !contact condition
                end Do !IP2
             end DO
          end DO
       end DO
999    CONTINUE
    end DO

!!$    PRINT*,'MAXIMUM OVERLAP = ', MAX_OVERLAP
!!$    IF(L_MAX.NE.0)THEN 
!!$       DIST(:) = DES_POS_NEW(L_MAX, :) - DES_POS_NEW(PNO_MAX,:)
!!$       DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
!!$
!!$       PRINT*,  L_MAX, PNO_MAX, DISTMAG, DES_RADIUS(L_MAX)+DES_RADIUS(PNO_MAX)
!!$       PRINT*,'MAXIMUM OVERLAP, PART POSs L', DES_POS_NEW(L_MAX,:)
!!$       PRINT*,'MAXIMUM OVERLAP, PART POSs J', DES_POS_NEW(PNO_MAX,:)
!!$       PRINT*,'MAXIMUM OVERLAP, PART CELLS L', PIJK(L_MAX,:)
!!$       PRINT*,'MAXIMUM OVERLAP, PART CELLS J', PIJK(PNO_MAX,:)
!!$    end IF

!    goto 100 


    PRINT*,'NOW DELETING PARTICLES CLOSER THAN ', MIN_PART_SEP
    DO LL = 1, PARTICLES
       IF(IFPON(LL).EQ.0) GOTO 9991
       
       PC(:) =  PIJK(LL,1:3)
       II = PC(1)
       JJ = PC(2)
       KK = PC(3)
       IP1 = II
       IM1 = II
       JP1 = JJ
       JM1 = JJ
       KM1 = KK
       KP1 = KK
!!$
       IP1 =II+1! MIN(IMAX1, II+1)
       IM1 =II-1! MAX(IMIN1, II-1)
       JP1 =JJ+1! MIN(JMAX1, JJ+1)
       JM1 =JJ-1! MAX(JMIN1, JJ-1)
         
       IF(DIMN.EQ.3) THEN 
          KP1 =KK+1! MIN(KMAX1, KK+1)
          KM1 =KK-1! MAX(KMIN1, KK-1)
       end IF
       !PRINT*,
       !PRINT*,'PC = ', IP1, JP1, KP1

       DO KK = KM1, KP1
          DO JJ = JM1, JP1
             DO II = IM1, IP1
                
                I = II
                J = JJ
                K = KK

                XPER_FAC = 0
                YPER_FAC = 0
                ZPER_FAC = 0
                PER_COND = .FALSE.
                IF(II.GT.IMAX1) THEN 
                   IF(INTX_PER) THEN 
                      I = IMIN1
                      XPER_FAC = one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true EAST',I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      I = IMAX1
                   ENDIF
                ENDIF

                IF(II.LT.IMIN1) THEN 
                   IF(INTX_PER) THEN 
                      I = IMAX1
                      XPER_FAC = -one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true WEST', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE 
                      I = IMIN1
                   ENDIF
                ENDIF

                IF(JJ.GT.JMAX1) THEN 
                   IF(INTY_PER) THEN 
                      J = JMIN1
                      YPER_FAC = one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true NORTH', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      J = JMAX1
                   ENDIF
                ENDIF

                IF(JJ.LT.JMIN1) THEN 
                   IF(INTY_PER) THEN 
                      J = JMAX1
                      YPER_FAC = -one
                      PER_COND = .true.
                      !WRITE(*,*) 'cond true SOUTH', I,J,K,SIZE(PIC(I,J,K)%p)
                   ELSE
                      J = JMIN1
                   ENDIF
                ENDIF

                IF(DIMN.EQ.3) THEN 
                   IF(KK.GT.KMAX1) THEN 
                      IF(INTZ_PER) THEN
                         K = KMIN1
                         ZPER_FAC = one
                      ELSE
                         K = KMAX1
                      ENDIF
                   ENDIF

                   IF(KK.LT.KMIN1) THEN 
                      IF(INTZ_PER) THEN
                         K = KMAX1
                         ZPER_FAC = -one
                      ELSE
                         K = KMIN1
                      ENDIF
                   ENDIF
                ENDIF

                If (ASSOCIATED(PIC(I,J,K)%p)) then
                   NPG = SIZE(PIC(I,J,K)%p)
                Else
                   NPG = 0
                Endif
!!$                IF(PER_COND) THEN 
!!$                   WRITE(*,*) 'pC = ', pc
!!$                   WRITE(*,*) 'II, JJ = ', II, JJ
!!$                   WRITE(*,*) 'I, J = ', I, J
!!$                   WRITE(*,*) 'XYPER_FAC ', XPER_FAC, YPER_FAC
!!$                   WRITE(*,*) 'DES_VEL_NEW = ', DES_POS_NEW(PNO,:)
!!$                   WRITE(*,*) 'MODIFIED POSITION = ', CORD_PNO(:)
!!$                   !READ(*,*)
!!$                ENDIF

                
                Do IP2 = 1,NPG
                   PNO = PIC(I,J,K)%p(ip2)
                   IF(PNO.GT.LL.AND.IFPON(PNO).EQ.1) THEN 
                         



                         
                      !  PRINT*,'LL, PNO = ', LL, PNO, NPG, I,J,K
                      R_LM = DES_RADIUS(LL) + DES_RADIUS(PNO)
                      !R_LM = FACTOR_RLM*R_LM
                      R_LM = R_LM + MIN_PART_SEP
                      CORD_PNO(1) = DES_POS_NEW(PNO,1) + XPER_FAC*(LX)
                      CORD_PNO(2) = DES_POS_NEW(PNO,2) + YPER_FAC*(LY+ONE)
                      IF(DIMN.EQ.3) THEN 
                         CORD_PNO(3) = DES_POS_NEW(PNO,3) + ZPER_FAC*(LZ+ONE)
                      ENDIF

                      DIST(:) = CORD_PNO(:) - DES_POS_NEW(LL,:)
                      DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
                      
                      !PRINT*,'D =', DISTMAG, R_LM
                      ALREADY_NEIGHBOURS = .FALSE.
                      
                      !IF(LL.EQ.1.AND.PNO.EQ.2) THEN 
                   !   WRITE(*,*)'CORD=', CORD_PNO(1)
                      !   WRITE(*,*)'DISTMAG', DISTMAG, R_LM-DISTMAG
                      !ENDIF
                      
                      
                      IF( DISTMAG.LT.R_LM) THEN 
                         !IF(PER_COND)THEN 
!!$                         
                            !PRINT*,'PARTICLE DELETED DUE TO MIN_PART_SEP', PER_COND
                            !PRINT*,'LL, PNO = ', LL, PNO
                            !PRINT*,'POS = ', DES_POS_NEW(LL,:), DES_POS_NEW(PNO,:)                         
                         !end IF
                         
                         ifpon(PNO) = 0
                      end IF
                      
                   end IF !contact condition
                end Do !IP2
             end DO
          end DO
       end DO
9991   CONTINUE
    end DO

    !STOP

100 continue
    
    

    DEALLOCATE(PIJK, PIC, CND, NEIGHBOURS)


  END SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH_MAT

  REAL(prcn) FUNCTION epan(x,h)
    IMPLICIT NONE
    REAL(prcn):: x, h
    epan=0.0
    IF ( (x .GE. (-h)) .AND. (x .LE. h) )THEN
       epan = (0.75/h)*(1.-x*x/(h*h))
    ENDIF
    RETURN
  END FUNCTION epan


END MODULE maternmod
    
 
