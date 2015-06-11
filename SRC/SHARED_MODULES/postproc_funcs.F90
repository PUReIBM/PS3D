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

MODULE postproc_funcs
  USE precision
  USE constants
  USE global_data, ONLY : ndim, RUN_NAME, lybyd, dbydx, dia_phys, minunitno, maxunitno
  use general_funcs
  Implicit NONE 
  Type ::  vertindex
     Private
     integer, dimension(3):: vindex
  end Type vertindex
  
  Type ::  cellindex
     Private
     integer, dimension(3):: cindex
  end Type cellindex

  Type(vertindex), Dimension(:), Allocatable, PRIVATE  :: gridvertno
  Type(cellindex), Dimension(:), Allocatable, PRIVATE  :: gridcellno
  
  integer(8) :: NVERT
	integer :: NCELL

CONTAINS 

  SUBROUTINE calc_gofr(n, xc, rad, my, mbox, xperiodic, nrbins, ndim, rescaling, gofr, rho_est, radbin)

    USE precision
    USE constants 
    USE general_funcs
    implicit none

    integer, Intent(in) :: ndim,mbox,my,n, nrbins
    real*8, Intent(in) ::  xc(n,ndim), rad(n)
    real*8, Intent(out) ::  gofr(nrbins), rho_est(nrbins), radbin(nrbins)

    Logical, INTENT(in) ::  xperiodic, rescaling
    real*8 epan,rmax_part,c

    external epan
    REAL*8 rmin,rmax,dr, xp(n,ndim)
    real*8 h,numdens_M,dij,Adenom,W1min,rmax_contrib
    real*8 rmin_contrib, L(ndim),Lint(ndim),r,temp


    integer ir,i,idim,j,irmin, irmax, ifpon(n), gofunit, np

    xp(1:n, 1:ndim) = xc(1:n, 1:ndim) 


    !open(gofunit,file='gofr.dat',status='unknown')
    !	endif
    !write(gofunit,*)'Zone'

    c = 0.15d0

    rmin = zero
    rmax = zero 
    DO idim = 1, ndim
       rmax = rmax + one!my**2
    ENDDO
    !rmax = SQRT(three)  
    rmax = SQRT(rmax)
    dr = (rmax - rmin)/float(nrbins-1)

    do ir = 1,nrbins
       rho_est(ir) = 0.0
    enddo
    rmax_part = MAXVAL(rad(1:n))

    do i=1,n
       do idim=1,ndim

          IF(rescaling) THEN
             if(idim.eq.1) then
                if(.not.xperiodic)then
                   xp(i,idim) = xp(i,idim)-rmax_part-3
                endif
                xp(i,idim) = (xp(i,idim)-1)/(mbox-1)
             else 
                xp(i,idim) = (xp(i,idim)-1)/(my)
             endif
          end IF !rescaling 

       enddo
       ifpon(i) = 1
    enddo

    do idim=1,ndim
       L(idim)= one!my
    enddo
    
    np = n
    numdens_M = n/(L(1)*L(2)*L(3))             !since L is unity

    !calculate bandwidth

    h = c/sqrt(numdens_M)

    write(*,*)'Bandwidth...',h

    !transfer positions of particles to xp matrix


    !Start calculation of rho^2(r) and g(r)=rho^2(r)/lamda^2
    !separate calculation of rho_2

    do i = 1, np

       if(mod(i,100).eq.0)write(*,*)i,' particles done..'

       do j = 1, np
          if ( j.ne.i .and.(ifpon(i).eq.1) .and. (ifpon(j).eq.1) ) then

             !calculate the inter-point distance dij

             dij = 0.0
             Adenom = 1.0
             do idim = 1, ndim
                dij = dij + (xp(i,idim) - xp(j,idim))**2

                !calculate contributions to the denominator A(Wxj \int Wxi)
                !our window is a unit cube

                !calculate the intersection of the segments (xp(i,idim), xp(i,idim)
                !+ L(idim)) and (xp(j,idim), xp(j,idim)+L(idim))

                !Find global min of Wi and Wj (in x,y,z succ), in x call this xmin;
                !if xp(i, ) = xmin we are in i, else in j; call this window 1
                !find the diff betn max window 1 - min window 2
                !if this diff is positive, there is a nonzero intersection,
                !else this diff is negative and the windows don't intersect

                W1min = min(xp(i,idim),xp(j,idim))
                if ( xp(i,idim) .eq. W1min ) then
                   !the window 1 is the i window
                   !Lint is the length of the intersecting segment in dim idim
                   !L(idim) = 1.0 in our case
                   Lint(idim) = max(0.d0,(xp(i,idim)+L(idim)-xp(j,idim)))
                else
                   !the window 1 is the j window
                   Lint(idim) = max(0.d0,(xp(j,idim)+L(idim)-xp(i,idim)))
                endif
                Adenom = Adenom*Lint(idim) !find intersecting area
                if(Adenom.eq.0.0)then
                   write(*,*)'Adenom...equal to zero..',i,j
                   stop
                endif
             enddo
             dij = sqrt(dij)

             !Use the Epanecnikov kernel; note that each interpoint distance
             !contributes to r bins that are within the kernel width +/- h

             !kernel width parameter h is recommended as c/sqrt(lambda)
             !take c = 0.15

             !calculate rmin and rmax that this pair contributes to

             rmin_contrib = dij - h
             rmax_contrib = dij + h

             !for each bin r(ir) such that rmin <= r(ir) <=rmax calculate
             !the value of rho_est

             irmin = (rmin_contrib-rmin)/dr
             irmax = (rmax_contrib-rmin)/dr
             do ir = irmin, irmax
                r = rmin + ir*dr
                temp = r - dij
                rho_est(ir) = rho_est(ir) + epan(temp,h)/(Adenom)
             enddo
          endif
       enddo

    enddo

    do ir = 1, nrbins
       r = rmin + ir*dr
       if(ndim.eq.1)then
          rho_est(ir) = rho_est(ir)/(2.*r)
       elseif(ndim.eq.2)then
          rho_est(ir) = rho_est(ir)/(2.*pi*r)
       elseif(ndim.eq.3)then
          rho_est(ir) = rho_est(ir)/(4.*pi*r**2)
       endif
       gofr(ir)  = rho_est(ir)/(numdens_M**2)
       radbin(ir) = r 
       !write(gofunit,1000)r, rho_est(ir), rho_est(ir)/(numdens_M**2)
    enddo

    !1000 format(10(E20.10,1x))

    !close(gofunit, status="keep")
  end SUBROUTINE calc_gofr

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

subroutine calc_numdens(nbody,xc,radbdy,ndim, nbins,avg_meth, n_mean, n_var, vol_bin)
  
    USE precision
    USE constants 
    USE general_funcs
    implicit none
  
  integer, Intent(in) :: nbody, ndim
  integer :: i,ncount,j,n,nsim,avg_meth, nbins, ix, iy, iz, idim
  
  double precision :: asize,xi,yi,zi,rmin,rmax, dr 
  double precision,Intent(in) :: xc(nbody,ndim), radbdy(nbody)
  double precision, Intent(out) :: n_mean(nbins),n_var(nbins),vol_bin(nbins)
  double precision :: dist, cr_x, temp
  !double precision :: number_per_cell(nbinmax, nbinmax, nbinmax), dx
  !double precision :: vel_per_cell_x(nbinmax, nbinmax, nbinmax)
  !double precision :: vel_per_cell_y(nbinmax, nbinmax, nbinmax)
  !double precision :: vel_per_cell_z(nbinmax, nbinmax, nbinmax), dx


!!!define a cube with center at 0.5,0.5,0.5, with size 'a'
  
  !open(unit=77,file='sndord_traninv.dat',status='unknown')
  !open(unit=78,file='avg_num_vel.dat',status='unknown')
  
  n = nbody 
  
    rmin = MINVAL(radbdy(1:nbody))
    rmin = rmin/20.d0
    rmax = zero
    DO idim = 1, ndim
       rmax = rmax + one
    ENDDO
    rmax = sqrt(rmax)/two
    
    dr = (rmax-rmin)/nbins
    asize = zero

  do j=1,nbins
     ncount = 0
     asize = asize + dr    
     
     do i=1,n
        xi = xc(i,1)
        yi = xc(i,2)
        zi = xc(i,3)
        
        if(avg_meth.eq.1)then
	   dist = ((xi-0.5d0)**2 + (yi-0.5d0)**2 + (zi-0.5d0)**2)**0.5d0

	   if(dist.le.asize)then
	      ncount = ncount + 1
           endif
        
        elseif(avg_meth.eq.0)then

           cr_x  = 0.25

           if((xi.lt.asize+cr_x.and.xi.gt.cr_x-asize).and. &
               & (yi.lt.asize+0.5.and.yi.gt.0.5-asize).and. &
               & (zi.lt.asize+0.5.and.zi.gt.0.5-asize))then
               ncount = ncount + 1
           elseif((-(1.-xi).lt.asize+cr_x.and.-(1.-xi).gt.cr_x-asize).and. &
               & (yi.lt.asize+0.5.and.yi.gt.0.5-asize).and. &
               & (zi.lt.asize+0.5.and.zi.gt.0.5-asize))then
               ncount = ncount + 1
           end if
        
        end if

     end do
     
     if(avg_meth.eq.1)then
       n_mean(j) = dble(ncount)
       n_var(j) = dble(ncount)**2
       vol_bin(j) = 4.d0/three*pi*asize**3 
     endif
       ! write(77,'10(E20.10,1x)')(2.*asize)**3, dble(ncount), dble(n)*(2.*asize)**3., & 
        !  &  dble(ncount)**2, (dble(n)*(2.*asize)**3.)**2
     !elseif(iopt.eq.1)then
      !  write(77,'10(E20.10,1x)')4./3.*pi*asize**3, dble(ncount), dble(n)*4./3.*pi*asize**3., & 
       !   &  dble(ncount)**2, (dble(n)*4./3.*pi*asize**3.)**2
     !endif
     
  end do

!!!! compute the mean of n and u for each MIS

  !nbins = 15 

  !dx = 1./dble(nbins)

  !number_per_cell(1:nbins,1:nbins,1:nbins) = 0
  !vel_per_cell_x(1:nbins,1:nbins,1:nbins) = 0.d0
  !vel_per_cell_y(1:nbins,1:nbins,1:nbins) = 0.d0
  !vel_per_cell_z(1:nbins,1:nbins,1:nbins) = 0.d0

  !do i=1,n
   !  ix = int(x(i)/dx) + 1
    ! iy = int(y(i)/dx) + 1
    ! iz = int(z(i)/dx) + 1

     !number_per_cell(ix,iy,iz) = number_per_cell(ix,iy,iz) + 1

     !vel_per_cell_x(ix,iy,iz) = vel_per_cell_x(ix,iy,iz) + vx(i)
     !vel_per_cell_y(ix,iy,iz) = vel_per_cell_y(ix,iy,iz) + vy(i)
     !vel_per_cell_z(ix,iy,iz) = vel_per_cell_z(ix,iy,iz) + vz(i)

  !end do

 ! write(78,*)'Zone i = ',nbins,' j = ',nbins,' k = ',nbins

  !do iz = 1, nbins
   !  do iy = 1, nbins
    !    do ix = 1, nbins
!	   if(number_per_cell(ix,iy,iz).eq.0)then
!	      temp = 1.
 !          else
  !            temp = 1./dble(number_per_cell(ix,iy,iz))
!	   endif
 !          write(78,'10(E20.10,1x)')dble(ix),dble(iy),dble(iz),dble(number_per_cell(ix,iy,iz)), &
  !              & dble(number_per_cell(ix,iy,iz))/dx**3., & 
   !             & vel_per_cell_x(ix,iy,iz)*temp, vel_per_cell_y(ix,iy,iz)*temp, vel_per_cell_z(ix,iy,iz)*temp
!
 !          write(91,'15(E20.10,1x)')dble((ix-1) + (iy-1)*nbins + (iz-1)*nbins*nbins+ 1), &
  !               & dble(number_per_cell(ix,iy,iz)), &
   !             & dble(number_per_cell(ix,iy,iz))/dx**3., &
    !            & vel_per_cell_x(ix,iy,iz)*temp, vel_per_cell_y(ix,iy,iz)*temp, vel_per_cell_z(ix,iy,iz)*temp
     !   end do
     !end do
  !end do

end subroutine calc_numdens

SUBROUTINE calc_vector_correlation(n, xc, rad, my, mbox, xperiodic, nrbins,&
     & ndim, rescaling, fij, radbin, vector)
  ! computes the correlation of a vector in a periodic box. This is a
  ! rank 2 tensor.
  USE precision
  USE constants
  IMPLICIT NONE
  integer, Intent(in) :: ndim,mbox,my,n, nrbins
  real*8, Intent(in) ::  xc(n,ndim), rad(n), vector(n,ndim)
  real*8, Intent(inout) ::  fij(ndim,ndim,nrbins), radbin(nrbins)
  
  Logical, INTENT(in) ::  xperiodic, rescaling
  real*8 rmax_part
  integer gofr(nrbins)
  REAL*8 rmin,rmax,dr, xp(n,ndim), tempr(ndim)
  real*8 :: L(ndim),Lint(ndim),r,temp, vmean(ndim),vtilde(n,ndim),vvar(ndim)
  
  
  integer i,idim,j,ibin, dimi, dimj
  
  xp(1:n, 1:ndim) = xc(1:n, 1:ndim) 

  do ibin = 1,nrbins
     fij(:,:,ibin) = 0.0
     radbin(ibin) = 0.0
     gofr(ibin) = 0
  enddo

  do idim=1,ndim
     vmean(idim) = SUM(vector(1:n,idim))/n
     vtilde(:,idim) = vector(:,idim)-vmean(idim)
     vvar(idim) = zero
  end do
  
  
  rmax_part = MAXVAL(rad(1:n))
  
  do i=1,n
     do idim=1,ndim
        
        IF(rescaling) THEN
           if(idim.eq.1) then
              if(.not.xperiodic)then
                 xp(i,idim) = xp(i,idim)-rmax_part-3
              endif
              xp(i,idim) = (xp(i,idim)-1)/(mbox-1)
           else 
              xp(i,idim) = (xp(i,idim)-1)/(my)
           endif
        end IF !rescaling 
        
     enddo

  enddo

  do idim=1,ndim
     L(idim)= one!my
  enddo

  rmin = zero
  rmax = zero 

!!$  rmax = DSQRT(3.d0)*(L(1))!/two 
  rmax = L(1)/two
  dr = (rmax - rmin)/float(nrbins-1)  

  do i = 1,n-1
     do j = i+1, n
        r = zero ! separation
        do idim = 1, ndim
           tempr(idim) = xp(i,idim) - xp(j,idim) ! compute the separation
           ! in each dimension
           if((ABS(tempr(idim))).gt.rmax) then
              if(tempr(idim).lt.zero) then
                 tempr(idim) = tempr(idim) + L(idim)
              else
                 tempr(idim) = tempr(idim) - L(idim)
              end if
           end if
           r = r + tempr(idim)**2.d0
        end do
        r = DSQRT(r)
        if(r.gt.rmax) goto 1001
        
        ibin = (r-rmin)/dr + 1
!!$        PRINT*, 'ibin = ', ibin
        gofr(ibin) = gofr(ibin)+1
        do dimj = 1, ndim
           do dimi = 1, ndim
              fij(dimi,dimj,ibin) =  fij(dimi,dimj,ibin) + vtilde(i&
                   &,dimi)*vtilde(j,dimj)
           end do
        end do
        
1001    continue
        
     end do 
     vvar(:) = vvar(:) + vtilde(i,:)*vtilde(i,:)
  end do
  vvar(:) = vvar(:) + vtilde(n,:)*vtilde(n,:)
  vvar(:) = vvar(:)/real(n,prcn)
  
  do ibin=1,nrbins
     radbin(ibin) = rmin + dr*(ibin-1)
     do dimj = 1, ndim
        do dimi = 1, ndim
           if(gofr(ibin).gt.0)then
              fij(dimi,dimj,ibin) = (fij(dimi,dimj,ibin)/gofr(ibin))/vvar(dimi)
           end if
           
        end do
     end do
     
  end do
!!$  open(unit=2002,file='gofr.dat',form='formatted',status='unknown')
!!$  
!!$  do ibin=1,nrbins
!!$     write(2002,'((2x,f12.8),2x,i3)') radbin(ibin), gofr(ibin)
!!$  end do
!!$  close(2002,status='keep')
  
end SUBROUTINE calc_vector_correlation

#if 0
SUBROUTINE scalar_two_point_correlation_pencil(phir,fij,radbin, conf_l95, conf_r95, nrbins)
  USE precision
  USE constants  
  !USE global_data
  !USE scalar_data
  IMPLICIT NONE
  integer :: i,j,k,i2,j2,k2, nrbins, fluid_pts, isp, ibin, idim, fluidpts_plane(mx1), jk_count
  real(prcn) , Intent(out), dimension(:) :: fij, radbin
  
  real(prcn) , Intent(in), dimension(:,:,:,:) :: phir
  real(prcn) :: rmax, rmin, drbin, tempr, L(3), vvar(nspmx), vtilde(2,nspmx), xp(2,3), gofr(nrbins),fij_loc(my*mz,nrbins), phi_fluid_mean_plane(mx1,nspmx), phimodmean_plane(mx1,nspmx),conf_l95(nrbins),conf_r95(nrbins),conf_l99(nrbins),conf_r99(nrbins)
  real(prcn) :: vvar_plane(my,mz,nspmx)
  

  
  fluid_pts = 0
  L(:) = doml(:)
  rmax = L(1)/two
  rmin = zero 
  drbin = (rmax - rmin)/float(nrbins-1)  
  PRINT*,'L =', L(:)
  WRITE(*,'(A25,3(2x,g12.5))')'rmin, rmax, drbin', rmin, rmax, drbin
  phi_fluid_mean = zero 
  phi_solid_mean = zero 
  phimodmean = zero 
  phi_fluid_mean_plane = zero 
  phimodmean_plane = zero 
  fluidpts_plane = 0
  do ibin = 1,nrbins
     fij(ibin) = 0.0
     radbin(ibin) = 0.0
     gofr(ibin) = 0
     fij_loc(:,ibin) = zero
     conf_l95(ibin) = zero 
     conf_r95(ibin) = zero
     conf_l99(ibin) = zero
     conf_r99(ibin) =  zero
  enddo
  
  
  vvar = zero 
  isp = 1
  DO k = 1, mz 
     DO j = 1, my
        DO i = 1, mx
           if(i.le.mx1) then 
              if(fluid_atijk(i,j,k)) then 
                 phi_fluid_mean(isp) = phi_fluid_mean(isp) + phir(i,j,k,isp)
                 phi_fluid_mean_plane(i,isp) = phi_fluid_mean_plane(i,isp) + phir(i,j,k,isp)
                 fluidpts_plane(i) = fluidpts_plane(i)+1
                 
              ELSE
                 phi_solid_mean(isp) = phi_solid_mean(isp)+ phir(i,j,k,isp)
              end if
              phimodmean(isp) = phimodmean(isp) +  phir(i,j,k,isp)
              phimodmean_plane(i,isp) = phimodmean_plane(i,isp) +  phir(i,j,k,isp)
              
           endif
     
        end DO
     end DO
  end DO
  phi_fluid_mean(1:nspmx) = phi_fluid_mean(1:nspmx)/(count_fluid)
  phi_solid_mean(1:nspmx) = phi_solid_mean(1:nspmx)/(count_solid)
  phimodmean = phimodmean/(mx1*my*mz)

  do i = 1, mx1 
     phi_fluid_mean_plane(i,1:nspmx) = phi_fluid_mean_plane(i,1:nspmx)/real(fluidpts_plane(i),prcn)
     
     phimodmean_plane(i,:) = phimodmean_plane(i,:)/real(my*mz,prcn)
     
     
  end do

     
     do k = 1, mz
        do j = 1, my 
           
           fluid_pts = 0
           vvar_plane(j,k,isp) = zero 
           do i = 1, mx1 
              if(fluid_atijk(i,j,k)) then 
                 fluid_pts = fluid_pts+1
                 vtilde(1,isp) = phir(i,j,k,isp)-phi_fluid_mean(isp)
                 
                 vvar_plane(j,k,isp) = vvar_plane(j,k,isp)+ vtilde(1,isp)*vtilde(1,isp)
              end if
              
           end do
           if(fluid_pts.gt.0) vvar_plane(j,k,isp) = vvar_plane(j,k,isp)/real(fluid_pts,prcn)
        end do
     end do
  
  
  WRITE(*,'(A26,3(2x,g12.5))') "PHI MEAN F and SF = ",  phi_fluid_mean(1:nspmx), phi_fluid_mean(1:nspmx)*(one-maxvolfrac)
  WRITE(*,'(A26,3(2x,g12.5))') "PHI MEAN S and SF = ",  phi_solid_mean(1:nspmx), phi_solid_mean(1:nspmx)*maxvolfrac  


  jk_count = 0
  i = 1
  DO k = 1, mz 
     DO j = 1, my
        gofr = 0
        jk_count = jk_count+1
        fij_loc(jk_count,:) = zero 
        DO i = 1, mx1
           ! i = 1
           if (.NOT.(fluid_atijk(i,j,k))) goto 2001
           
           xp(1,1) = i*dx
           xp(1,2) = j*dy
           xp(1,3) = k*dz
           !SHUD HAVE i-1, j-1. but i-i2 will still remain the same. 
           vtilde(1,isp) = phir(i,j,k,isp)-phi_fluid_mean(isp)
           j2 = j
           k2 = k 
           vvar = zero 
           fluid_pts = 0
           Do i2 = i, mx1,1
              
              if (.NOT.(fluid_atijk(i2,j2,k2))) goto 3001
              xp(2,1) = i2*dx
              xp(2,2) = j2*dy
              xp(2,3) = k2*dz
              
              vtilde(2,isp) = phir(i2,j2,k2,isp)-phi_fluid_mean(isp)
              
              tempr = xp(1,1)-xp(2,1)
              ! in each dimension
              if((ABS(tempr)).gt.rmax) then
                 if(tempr.lt.zero) then
                    tempr = tempr + L(1)
                 else
                    tempr = tempr - L(1)
                 end if
              end if
              
              r = abs(tempr)
              ibin =NINT(((r-rmin)/drbin))+1
              
              
              !if(j.eq.1..and.k.eq.1) WRITE(*,'(A,2(2x,i2),5(2x,g12.5))') 'i,i2,r,ibin',i,i2,r,r/drbin+1,((r-rmin)/drbin)+1, ibin

              if(ibin.le.0.or.ibin.gt.nrbins) PRINT*, 'ibin = ', ibin
              
              gofr(ibin) = gofr(ibin)+1
              fij_loc(jk_count,ibin) =  fij_loc(jk_count,ibin) + vtilde(1,isp)*vtilde(2,isp)
              
              !vvar(isp) = vvar(isp) + vtilde(2,isp)*vtilde(2,isp)

              
              fluid_pts = fluid_pts+1
              
3001          continue 
              
           end DO !I2
           
!           exit 
2001       continue 
           
        end DO ! I
        !PRINT*,'I = ', I
        !vvar(isp) = vvar(isp)/real(fluid_pts,prcn)
        do ibin=1,nrbins
           radbin(ibin) = rmin + drbin*(ibin-1)
! if(j.eq.1..and.k.eq.1) write(*,'(A,5(2x,g12.5))')'fij =', fij_loc(jk_count,ibin), vvar_plane(j,k,isp), gofr(ibin),(fij_loc(jk_count,ibin)/gofr(ibin))/vvar_plane(j,k,isp)
              if(gofr(ibin).gt.0)then
                 fij_loc(jk_count,ibin) = (fij_loc(jk_count,ibin)/gofr(ibin))/vvar_plane(j,k,isp)
              end if
              
        end do

     
     end do ! J
  end DO ! K
  
  !PRINT*,'COUNT_FLUID, FLUID_PTS = ', COUNT_FLUID, FLUID_PTS
  PRINT*,'JK_COUNT = ', JK_COUNT
  !jk_count = 1
  do ibin = 1, nrbins
     CALL calc_confidence(fij_loc(1:jk_count,ibin),jk_count,conf_l95(ibin),conf_r95(ibin),conf_l99(ibin),conf_r99(ibin),fij(ibin))
  end do
END SUBROUTINE scalar_two_point_correlation_pencil
#endif

SUBROUTINE calc_confidence(qstat,nreal,climl95,climr95,climl99,climr99,ave)
  USE precision
  Use constants
  IMPLICIT NONE
  REAL(prcn), DIMENSION(nreal), INTENT(in) :: qstat
  INTEGER, INTENT(in) :: nreal
  INTEGER :: jreal,istat,j
  REAL(prcn) , INTENT(out) :: ave, climr95,climl95,climl99,climr99
  REAL(prcn) :: ds, dave,dvar,ctemp,ratio,svar,sdev
  
  ds=0.d0
  DO jreal=1,nreal
     ds = ds+DBLE(qstat(jreal))
  ENDDO
  dave=ds/DBLE(float(nreal))
  ave = (dave)
  dvar=0.d0
  DO j=1,nreal
     ds=DBLE(qstat(j))-dave
     dvar=dvar+ds*ds
  ENDDO
  dvar=dvar/DBLE(float(nreal-1))
  svar=sngl(dvar)
  sdev=SQRT(svar)
  ctemp = sdev/SQRT(float(nreal))
  climl95 = ave - 2.492*ctemp
  climl99 = ave - 3.467*ctemp
  climr95 = ave + 2.492*ctemp
  climr99 = ave + 3.467*ctemp
  ratio = (climr95-climl95)*100./ave
END SUBROUTINE calc_confidence


#if 0
SUBROUTINE scalar_two_point_correlation(phir,fij,radbin, nrbins)
  USE precision
  USE constants  
  !USE global_data
  !USE scalar_data
  IMPLICIT NONE
  integer :: i,j,k,ijk, IJK2,i2,j2,k2, nrbins, fluid_pts, isp, ibin, idim
  real(prcn) , Intent(out), dimension(:) :: fij, radbin
  
  real(prcn) , Intent(in), dimension(:,:,:,:) :: phir
  real(prcn) :: rmax, rmin, drbin, tempr(3), L(3), vvar(nspmx), vtilde(2,nspmx), xp(2,3), gofr(nrbins)

  WRITE(*,*) 'IN SCALAR TWO PT CORRELATION'
  CALL initialize_gridvertindex(mx1, my, mz)

  
  fluid_pts = 0
  L(:) = doml(:)
  rmax = L(1)/two
  rmin = zero 
  drbin = (rmax - rmin)/float(nrbins-1)  
  PRINT*,'L =', L(:)
  PRINT*,'rmin, rmax, drbin', rmin, rmax, drbin
  phi_fluid_mean = zero 
  phi_solid_mean = zero 
  phimodmean = zero 
  isp = 1

  do ibin = 1,nrbins
     fij(ibin) = 0.0
     radbin(ibin) = 0.0
     gofr(ibin) = 0
     
  enddo
  
  vvar = zero 
  
  DO IJK = 1, NVERT
    
     CALL ind1t3(ijk,i,j,k)

     if(i.le.mx1) then 
        if(fluid_atijk(i,j,k)) then 
           phi_fluid_mean(isp) = phi_fluid_mean(isp) + phir(i,j,k,isp)
        ELSE
           phi_solid_mean(isp) = phi_solid_mean(isp)+ phir(i,j,k,isp)
        end if
        
        phimodmean(isp) = phimodmean(isp) +  phir(i,j,k,isp)
     endif
     
  end DO
  
  phi_fluid_mean(1:nspmx) = phi_fluid_mean(1:nspmx)/(count_fluid)
  phi_solid_mean(1:nspmx) = phi_solid_mean(1:nspmx)/(count_solid)
  phimodmean = phimodmean/(mx1*my*mz)
  WRITE(*,'(A26,3(2x,g12.5))') "PHI MEAN F and SF = ",  phi_fluid_mean(1:nspmx), phi_fluid_mean(1:nspmx)*(one-maxvolfrac)
  WRITE(*,'(A26,3(2x,g12.5))') "PHI MEAN S and SF = ",  phi_solid_mean(1:nspmx), phi_solid_mean(1:nspmx)*maxvolfrac  
  
  DO IJK = 1, NVERT-1
     PRINT*,'IJK = ', IJK, NVERT
     CALL ind1t3(ijk,i,j,k)
     if (.NOT.(fluid_atijk(i,j,k))) goto 2000

     xp(1,1) = i*dx
     xp(1,2) = j*dy
     xp(1,3) = k*dz
     !SHUD HAVE i-1, j-1. but i-i2 will still remain the same. 
     vtilde(1,isp) = phir(i,j,k,isp)-phi_fluid_mean(isp)

     
     DO IJK2 = IJK, NVERT
        CALL ind1t3(ijk2,i2,j2,k2)
        if (.NOT.(fluid_atijk(i2,j2,k2))) goto 3000
        
        xp(2,1) = i2*dx
        xp(2,2) = j2*dy
        xp(2,3) = k2*dz
        
        
        vtilde(2,isp) = phir(i2,j2,k2,isp)-phi_fluid_mean(isp)
        r = zero 
        do idim = 1, ndim
           tempr(idim) = (xp(1,idim)-xp(2,idim))
           ! in each dimension
           if((ABS(tempr(idim))).gt.rmax) then
              if(tempr(idim).lt.zero) then
                 tempr(idim) = tempr(idim) + L(idim)
              else
                 tempr(idim) = tempr(idim) - L(idim)
              end if
           end if
           r = r + tempr(idim)**2.d0
        end do
        r = DSQRT(r)
        if(r.gt.rmax) goto 1001
        
        ibin = (r-rmin)/drbin + 1
        
        !PRINT*, 'ibin = ', ibin
        gofr(ibin) = gofr(ibin)+1
        fij(ibin) =  fij(ibin) + vtilde(1,isp)*vtilde(2,isp)
        
1001    continue

3000 continue 

     end DO !IJK2
     fluid_pts = fluid_pts+1
     vvar(isp) = vvar(isp) + vtilde(1,isp)*vtilde(1,isp)
2000 continue 
     
  end do !IJK
  
  CALL ind1t3(nvert,i,j,k)
  if(fluid_atijk(i,j,k)) then 
     vtilde(1,isp) = phir(i,j,k,isp)-phi_fluid_mean(isp)
     
     vvar(isp) = vvar(isp) + vtilde(1,isp)*vtilde(1,isp)
     
     fluid_pts = fluid_pts+1
  end if
  
  PRINT*,'COUNT_FLUID, FLUID_PTS = ', COUNT_FLUID, FLUID_PTS
  vvar(isp) = vvar(isp)/real(count_fluid,prcn)
     
     
  do ibin=1,nrbins
     radbin(ibin) = rmin + drbin*(ibin-1)
     if(gofr(ibin).gt.0)then
        fij(ibin) = (fij(ibin)/gofr(ibin))/vvar(isp)
     end if
  
  end do
end SUBROUTINE scalar_two_point_correlation
#endif


	Subroutine initialize_gridvertindex(ni,nj,nk)
		Implicit None 
		INTEGER, INTENT(IN) :: ni, nj, nk
		integer :: i,j,k,ii
		allocate(gridvertno(ni*nj*nk))
		NVERT = ni*nj*nk
		do k = 1, nk
			do j = 1, nj
				do i = 1, ni
					ii = i + (j-1)*ni + (k-1)*ni*nj
					gridvertno(ii)%vindex(1) = i
					gridvertno(ii)%vindex(2) = j
					gridvertno(ii)%vindex(3) = k
				enddo
			enddo
		enddo
	End Subroutine initialize_gridvertindex

    Subroutine initialize_gridcellindex(ni, nj, nk)
      
    Implicit None 
    
    INTEGER, INTENT(IN) :: ni, nj, nk
    integer :: i,j,k,ii
    allocate(gridcellno((ni-1)*(nj-1)*(nk-1)))
    
    NCELL = (ni-1)*(nj-1)*(nk-1)
    do k = 1, nk-1
       do j = 1,nj-1
          do i = 1,ni-1
             ii = (i-1)*(nj-1)*(nk-1) + (j-1)*(nk-1) + k 
             gridcellno(ii)%cindex(1) = i
             gridcellno(ii)%cindex(2) = j
             gridcellno(ii)%cindex(3) = k
          enddo
       enddo
    enddo
  End Subroutine initialize_gridcellindex

  Subroutine ind_cell1t3(ii,i,j,k)
    Implicit None 
    integer, intent(in) :: ii
    integer, intent(out) :: i,j,k
    i = gridcellno(ii)%cindex(1)
    j = gridcellno(ii)%cindex(2)
    k = gridcellno(ii)%cindex(3)
  End Subroutine ind_cell1t3

	Subroutine ind1t3(ii,i,j,k)
		Implicit None 
		integer(8), intent(in) :: ii
		integer, intent(out) :: i,j,k
		i = gridvertno(ii)%vindex(1)
		j = gridvertno(ii)%vindex(2)
		k = gridvertno(ii)%vindex(3)
	End Subroutine ind1t3

	subroutine vertindex_ijk(ii,i,j,k,nx,ny,nz)
		implicit none
		integer, intent(in) :: i, j, k, nx, ny, nz
		integer, intent(out) :: ii

		ii = i + (j-1)*nx + (k-1)*ny*nz
	end subroutine vertindex_ijk

#if 0
  SUBROUTINE calc_anisotropy(aij, zi, eta)
    IMPLICIT NONE

    REAL(prcn),Intent(in) :: aij(ndim,ndim)
    REAL(prcn),Intent(out) :: zi, eta
    Integer :: m,i,j,k
    Real(prcn) :: trace, bij(ndim,ndim), tmp_arr1(ndim,ndim),tmp_arr2(ndim,ndim),tmp,delta_ij

    trace = aij(1,1)+aij(2,2)+aij(3,3)
    !PRINT*,' trace =', trace 
    if(.not.(ABS(trace).GT.SMALL_NUMBER)) then
       PRINT*,'TRACE OF THE TENSOR IS VERY SMALL', trace
       RETURN
    end if
    do i = 1, ndim
       do j = 1, ndim
          if(i.eq.j)then 
             delta_ij = 1
          else
             delta_ij = 0
          endif
          bij(i,j) = aij(i,j)/(trace) - delta_ij/3.d0

       enddo

    enddo

    eta = zero
    zi = zero
    
    do i = 1,ndim
       do j = 1, ndim
          eta = eta + bij(I,J)*bij(j,i)
          do k= 1, ndim
             zi = zi + bij(i,j)*bij(j,k)*bij(k,i)
          end do
       end do
    end do

    !tmp_arr1 = matmul(bij,bij)
    
    !eta  = sqrt( (tmp_arr1(1,1) + tmp_arr1(2,2) + tmp_arr1(3,3))/6.d0 )
    eta = sqrt(eta/6.d0)
    zi = zi/6.d0
!    tmp_arr2 = matmul(bij, tmp_arr1) 
 !   zi = ( tmp_arr2(1,1) + tmp_arr2(2,2) + tmp_arr2(3,3) )/6.d0
    tmp = abs(zi)**(1./3.)
    zi = sign(tmp, zi)
    
  end subroutine calc_anisotropy
#endif


SUBROUTINE histogram(u,wt,n,nddu,bldd,brdd,hist)
  !nddu: number of bins for pdf formation
  !bldd: lefh hand side limit .... output
  !brdd: right side limie.... otuput
  USE precision
  USE constants
  USE randomno
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n , nddu
  REAL(prcn), INTENT(in), DIMENSION(:) :: u, wt
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: u_t
  
  REAL(prcn),  INTENT(out) :: bldd, brdd
  
  REAL(prcn),  INTENT(inout), DIMENSION(:) ::  hist(nddu)
    REAL(prcn) ::  xdiff
    
    REAL(prcn) :: vlmt, ave, adev, sdev, var, skew, curt
    
    INTEGER :: i, ibin
    bldd= 1.e25
    brdd=-1.e25
    ALLOCATE(u_t(n))
    CALL moment1(4, n, n, wt, u, ave,adev,sdev,var,skew&
         &,curt)
!!$    PRINT*,'ave,var,sdev, skew, curt=', ave, var,sdev, skew,curt
!!$    WRITE(*,*)'number of bins in hist..',nddu

  DO i=1,n
     u_t(i)=(u(i)-ave)/sdev
  ENDDO

  CALL moment1(4, n, n, wt, u_t, ave,adev,sdev,var,skew&
       &,curt)
!!$  PRINT*,'ave,var, sdev,skew, curt=', ave, var,sdev, skew,curt

    DO i=1,nddu
       hist(i)=0.0  
    ENDDO
    bldd = MIN(MINVAL(u_t(:)), bldd)
    brdd = MAX(MAXVAL(u_t(:)), brdd)
!!$    DO i=1,n
!!$       bldd = amin1(u(i),bldd)
!!$       brdd = amax1(u(i),brdd)
!!$    ENDDO

    xdiff  = (brdd - bldd)/float(nddu-1)

    DO i=1,n
       ibin = (u_t(i) - bldd)/xdiff + 1
       hist(ibin)=hist(ibin) + wt(i)/xdiff
    ENDDO

    DEALLOCATE(u_t)

  END SUBROUTINE histogram

SUBROUTINE plothist(hist,lb,ub,nbins,iunit,t,tref)
  USE precision
  USE constants
  IMPLICIT NONE
  REAL(prcn), INTENT(in), DIMENSION(:) ::hist
  
  INTEGER, INTENT(in) ::iunit,nbins
  REAL(prcn) ::  sum_x,lb,ub,dx, t, tref, tmp, tmp2
  INTEGER(8) :: i
  
  sum_x = lb
  dx= (ub-lb)/(float(nbins)-1)
  tmp = one/sqrt(twopi)
  WRITE(iunit,*)'Zone t="',t/tref,'"'
  DO i=1,nbins
     tmp2  = sum_x+dx/2
     WRITE(iunit,*)tmp2,hist(i)!, tmp*exp(-(tmp2*tmp2)/two)
     sum_x = sum_x + dx
  ENDDO
  
  RETURN
END SUBROUTINE plothist

SUBROUTINE compute_aivj(aivj,vj,ai,n)
  IMPLICIT NONE
  Integer, Intent(in) :: n
  REAL(prcn),Intent(inout) :: aivj(ndim,ndim)
  REAL(prcn),Intent(in) :: vj(n,ndim),ai(n,ndim)
  REAL(prcn) :: vivj(ndim,ndim)
  Integer :: m,i,j
  
  aivj = zero
  do i = 1, ndim
     do j = 1, ndim
        do m = 1, n
           aivj(i,j) = aivj(i,j) + ai(m,i)*vj(m,j)
        enddo
     enddo
  enddo
  aivj = aivj/real(n,prcn)
  !  vivj = vivj/n
  ! CALL calc_anisotropy(vivj)
end SUBROUTINE compute_aivj



	SUBROUTINE calculate_gofr_homog(n,xc, contact, grid, nrbins, rescaling, gofr, radbin, overlap)
		IMPLICIT NONE

		Integer, Intent(in) :: grid(ndim), n, nrbins
		Real(prcn), Intent(in) ::  xc(n,ndim)
		logical, intent(inout) :: contact(n,n)
		real(prcn), intent(inout) :: overlap
		real(prcn), Intent(out) ::  gofr(nrbins), radbin(nrbins)

		Logical, INTENT(in) ::  rescaling

		Integer :: i,idim,j,ibin, dimi, dimj, count

		REAL(prcn), ALLOCATABLE, DIMENSION(:,:) :: xp
		REAL(prcn) ::  L(ndim), volbin, rmax, rmin, dr, tempr(ndim), decimal, r

		ALLOCATE(xp(n,ndim))
    
		radbin(:) = zero
		gofr(:) = zero
    
!    if(rescaling)then
!       do i=1,n
!          do idim=1,ndim
!             if(idim.eq.1)then
!                xp(i,idim) = (xc(i,idim)-1)/(mbox-1)
!             else
!                xp(i,idim) = (xc(i,idim)-1)/(my)
!             end if
!          enddo
!       enddo
!    end if
    
!!$    DO idim=1,ndim
!!$       CALL uni_dist(xp(1:n,idim))
!!$    ENDDO !generate np uniformly 

!    do idim=1,ndim
!       L(idim)= one
!    enddo

		l(1) = dble(grid(1))
		l(2) = dble(grid(2))
		l(3) = dble(grid(3))

		rmin = zero
		rmax = L(1)/two
		dr = (rmax - rmin)/float(nrbins)

		count = 0
		do i = 1, n-1
			do j = i+1, n
				r = zero ! separation
				do idim = 1, ndim
					tempr(idim) = Abs(xc(i,idim) - xc(j,idim)) ! compute the separation

					if(tempr(idim).gt.rmax) then
						tempr(idim) = 2*rmax - tempr(idim)
					end if
					r = r + tempr(idim)**2.d0
				end do
				r = DSQRT(r)
				if(r.gt.rmax) goto 1001

				if  (r<=dia_phys) then
					if (overlap<dia_phys-r) overlap = dia_phys-r

					count = count+1
					contact(i,j) = .true.
					contact(j,i) = .true.
				endif

				ibin = (r-rmin)/dr + 1
				if(ibin>nrbins) ibin = nrbins

				gofr(ibin) = gofr(ibin) + one
1001			continue
			end do
		end do

		do ibin=1,nrbins
			radbin(ibin) = rmin + dr*(ibin-half)
			volbin = 4.d0*pi*(radbin(ibin)**2.d0)*dr
			gofr(ibin) = gofr(ibin)/real(n,prcn)
			gofr(ibin) = gofr(ibin)*L(1)*L(2)*L(3)*two/(n*volbin)
		end do
	end subroutine calculate_gofr_homog


	subroutine calc_gofr_2d(n, xc_2d, nx, ny, nbin, rad_2d, radbin, gr)
		implicit none
		integer, intent(in) :: n, nbin, nx, ny
		real(prcn), intent(in) :: xc_2d(n,2), rad_2d(n)
		real(prcn), intent(out) :: gr(nbin), radbin(nbin)

		integer :: i, j, ibin
		real(prcn) :: dist, dist1, dist2, dr, rmin, rmax, volbin

		rmin = zero
		rmax = ny/two
		dr = (rmax - rmin)/float(nbin)

		radbin(:) = zero
		gr(:)     = zero

		do i=1, n-1
			do j=i+1, n
				dist1 = abs(xc_2d(i,1)-xc_2d(j,1))
				if (dist1>rmax) dist1 = 2*rmax-dist1

				dist2 = abs(xc_2d(i,2)-xc_2d(j,2))
				if (dist2>rmax) dist2 = 2*rmax-dist2

				dist = sqrt(dist1**2 + dist2**2)
				if (dist>rmax) goto 100

				ibin = dist/dr+1

				if (ibin>nbin) ibin = nbin

				gr(ibin) = gr(ibin)+1
100			continue
			enddo
		enddo

		do ibin=1,nbin
			radbin(ibin) = rmin + dr*(ibin-half)
			volbin = 2*pi*radbin(ibin)*dr / (nx*ny)
			gr(ibin) = gr(ibin)/real(n,prcn)
			gr(ibin) = gr(ibin)/n / volbin *2/n
		end do
	end subroutine calc_gofr_2d

	SUBROUTINE GET_CONFIN(ndatain, confin)
		implicit none 
		INTEGER, INTENT(IN) :: ndatain 
		REAL(prcn), INTENT(OUT) ::  confin
		CHARACTER*50 :: table_name
		INTEGER :: tunit, ndata 
		LOGICAL :: filexist,isopen

		integer, parameter :: nmismax = 30, percmax = 11
		character*8 :: per_conf="95%"

		Type :: conf_interval
		 character*8 :: confper(percmax)
		 real(prcn) :: confi(percmax)
		END Type conf_interval

		TYPE(conf_interval) :: cis(nmismax)

		ndata = ndatain 
		!  WRITE(*,*)'IN GET CONFIN, NDATA = ', ndata
		table_name = "ci_lookup.dat"
		INQUIRE(FILE=trim(table_name),EXIST=filexist,OPENED=isopen)
		if(filexist)then
			tunit = getnewunit(minunitno,maxunitno)
			OPEN(unit=tunit,FILE=table_name,form='formatted')
			CALL initialize_lookup_table(tunit)
			if(ndata.gt.nmismax)ndata=nmismax
			if(ndata.gt.1)then
			confin = lookup_entry(ndata-1,per_conf)

			!WRITE(*,'(2(A,2x,g17.8))') 'DATA PTS = ', real(NDATA),' CONFIN = ', CONFIN
			else
			!WRITE(*,'(A,2x,g17.8)') 'ONLY ONE MIS: CONFIN = ', CONFIN
			confin = zero
			end if
			close(tunit,status='keep')

		else

			write (*,*) trim(table_name)//" DOES NOT EXIST. CHECK THE FILE...."
			stop
			confin = 1.812
		end if

	contains
		SUBROUTINE initialize_lookup_table(tunit)
			IMPLICIT NONE	


			CHARACTER*8 :: ciperc(percmax)
			INTEGER, Intent(in) :: tunit
			INTEGER :: iperc,isamples

			READ(tunit,*)(ciperc(iperc),iperc=1,percmax)

			do isamples=1,nmismax
				cis(isamples)%confper(1:percmax) = ciperc(1:percmax)
			end do
			do isamples=1,nmismax
				READ(tunit,*)(cis(isamples)%confi(iperc),iperc=1,percmax)
			end do

		END SUBROUTINE initialize_lookup_table

		RECURSIVE function lookup_entry(nsamples,entry) result(confis)
			IMPLICIT NONE
			INTEGER, Intent(in) :: nsamples
			CHARACTER*8, Intent(inout) :: entry
			REAL(prcn) :: confis
			Integer :: iperc

			do iperc = 1,percmax
				if(entry.eq.(cis(nsamples)%confper(iperc)))then
					confis = cis(nsamples)%confi(iperc)
					exit
				end if
			end do

			if(iperc.gt.percmax)then
			!       PRINT*,'ENTRY ', entry,' NOT FOUND. SETTING 95% CONFIDENCE INTERVAL'
				entry = TRIM('95%')
				confis = lookup_entry(nsamples,entry)
				RETURN
			end if
!    PRINT*,'CONFIDENCE INTERVAL FOR ', TRIM(entry),' CONFIDENCE AND ',nsamples, ' DEGREES OF FREEDOM IS ', confis 
			RETURN
		END function lookup_entry
	end SUBROUTINE GET_CONFIN


end MODULE postproc_funcs

real*8 function epan(x,h)
  implicit none
  real*8, Intent(in) :: x, h
  epan=0.0
  if ( (x .ge. (-h)) .and. (x .le. h) )then
     epan = (0.75d0/h)*(1.d0-x*x/(h*h))
  endif
  return
end function epan
      
      
