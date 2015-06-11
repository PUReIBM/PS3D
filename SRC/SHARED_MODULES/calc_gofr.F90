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
  Implicit NONE 
CONTAINS 

  SUBROUTINE calc_gofr(n, xc, rad, my, mbox, xperiodic, nrbins, ndim, rescaling, gofr, rho_est, rad_bin)

    USE PRECISION 
    USE constants 
    USE general_funcs
    implicit none

    integer, Intent(in) :: ndim,mbox,my,n, nrbins
    real*8, Intent(in) ::  xc(n,ndim), rad(n)
    real*8, Intent(out) ::  gofr(nrbins), rho_est(nrbins), rad_bin(nrbins)

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
                xp(i,idim) = (xp(i,idim)-1)/(my-1)
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
       rad_bin(ir) = r 
       !write(gofunit,1000)r, rho_est(ir), rho_est(ir)/(numdens_M**2)
    enddo

    !1000 format(10(E20.10,1x))

    !close(gofunit, status="keep")
  end SUBROUTINE calc_gofr

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

subroutine calc_numdens(nbody,xc,radbdy,ndim, nbins,avg_meth, n_mean, n_var, vol_bin)
  
  implicit none
  
  integer, Intent(in) :: nbody, ndim
  integer :: i,ncount,j,n,nsim,avg_meth, nbins, ix, iy, iz,
  
  double precision :: rsize,xi,yi,zi,rmin,rmax, dr 
  double precision,Intent(in) :: xc(nbody,ndim), radbdy(nbody)
  double precision, Intent(out) :: n_mean(nbins),n_var(nbins),vol_bin(nbins)
  double precision :: dist, cr_x, temp
  double precision :: number_per_cell(nbinmax, nbinmax, nbinmax), dx
  !double precision :: vel_per_cell_x(nbinmax, nbinmax, nbinmax)
  !double precision :: vel_per_cell_y(nbinmax, nbinmax, nbinmax)
  !double precision :: vel_per_cell_z(nbinmax, nbinmax, nbinmax), dx


!!!define a cube with center at 0.5,0.5,0.5, with size 'a'
  
  open(unit=77,file='sndord_traninv.dat',status='unknown')
  open(unit=78,file='avg_num_vel.dat',status='unknown')
  
  n = nbody 
  read(26,*)nsim
  
    rmin = MINVAL(radbody(1:nbody))/20.d0
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
     
     if(avg_meth.eq.0)then
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
end MODULE postporc_funcs

real*8 function epan(x,h)
  implicit none
  real*8, Intent(in) :: x, h
  epan=0.0
  if ( (x .ge. (-h)) .and. (x .le. h) )then
     epan = (0.75d0/h)*(1.d0-x*x/(h*h))
  endif
  return
end function epan
      
      
