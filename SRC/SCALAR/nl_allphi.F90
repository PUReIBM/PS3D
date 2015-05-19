MODULE nl_allphi
#include "../FLO/ibm.h"
  USE precision 
  USE constants 
  USE fftw_interface
  USE scalar_data 
  
  USE functions, ONLY : maxzero
#if !PARALLEL
  USE nlarrays, ONLY : uf1,uf2,uf3, phif1=>uf11, uf1phif=>uf12,&
       & uf2phif=>uf13,uf3phif=>uf23, ur1,ur2,ur3, phir1=>ur11,&
       & ur1phi=>ur12, ur2phi=>ur13, ur3phi=>ur23
#else
  USE nlarrays, ONLY : uf1,uf2,uf3, phif1=>uf11, uf1phif=>uf12,&
       & uf2phif=>uf13,uf3phif=>uf23, ur1,ur2,ur3, phir1=>ur11,&
       & ur1phi=>ur12, ur2phi=>ur13, ur3phi=>ur23, uatminus1, uatnxp2
#endif  
  IMPLICIT NONE 
  
  PRIVATE 
  INTEGER :: mx2
  PUBLIC :: form_nlphi
CONTAINS
  SUBROUTINE form_nlphi

    USE bcsetarrays, fdderiv=>fr  
    
    Use nlmainarrays, Only : ubcp, phirbcp=>nlbcp
    IMPLICIT NONE 
    INTEGER :: i,j,k, im1, ip1, isp,im2,ip2,umean_flag 
    REAL(prcn) :: Dplus_vj, Dminus_vj, max_speed,phifmeanloc(nspmx) &
         &,phismeanloc(nspmx), phimodmeanloc(nspmx), slope_factor, phi_ri
#if !PARALLEL
    REAL(PRCN) :: uplus(mx1), uminus(mx1), U_plus(mx1),&
         & U_minus(mx1), slope_u(mx1),slope_phi(mx1)
#else
    REAL(PRCN), ALLOCATABLE, DIMENSION(:) :: uplus, uminus, U_plus,&
         & U_minus, slope_u, slope_phi
#endif

#if 0
    if(irestart.eq.0) then 
       if(iglobstep.le.100) then
          slope_factor = nl_slope_factor*real(iglobstep-1,prcn)/100.d0
       ELSE
          slope_factor = nl_slope_factor
          !discr_scheme = "center"
       end if
    ELSE
       slope_factor = nl_slope_factor
    end if
#endif
    slope_factor = one
    phi_ri = one
    
    if(I_AM_NODE_ZERO)WRITE(*,'(A,2x,g17.8)')'SLOPE FACTOR IN NLPHI = ', slope_factor 
    
    mx2 = mx+1
    do isp = 1, nspmx
#if !PARALLEL
       do i = 1, nx !mx1
#else
       do i = 0, nx+1
#endif
          DO k = 1,mz
             DO j = 1, my2
!                if(i.eq.2) 
                IF(isp.eq.1.and.flow_converged) then
                   uf1(j,k) = u(i,j,k,1) 
                   
                   uf2(j,k) = u(i,j,k,2) 
                   
                   uf3(j,k) = u(i,j,k,3) 
                end IF
                
                phif1(j,k) = phif(i,j,k,isp)
             END DO
          END DO
          IF(isp.eq.1.and.flow_converged) then 
             CALL ff2cr(uf1,ubcp(i,:,:,1))
             CALL ff2cr(uf2,ubcp(i,:,:,2))
             CALL ff2cr(uf3,ubcp(i,:,:,3))
          end IF
          
          CALL ff2cr(phif1,phirbcp(i,:,:,isp))
          if(flow_converged) then
             umean_flag = 1
          else
             umean_flag = 0
          endif
          !IF flow is converged then umean has to be added to the transformed fluid velocity field, otherwise use the transformed in nl_flow earlier
          do k = 1, mz
             do j = 1, my 
                if(isp.eq.1) then
                   
                   IF(zero_flow_in_solid)  THEN 
                      
                      IF(Fluid_atijk(i,j,k)) THEN 
                         ubcp(i,j,k,1:3) = ubcp(i,j,k,1:3) + umean(1:3)*umean_flag
                      ELSE
                         ubcp(i,j,k,1:3) = zero
                      ENDIF
                   ELSE
                      ubcp(i,j,k,1:3) = ubcp(i,j,k,1:3) + umean(1:3)*umean_flag
                   end IF
                   
                end if
                
                phirbcp(i,j,k,isp) = phirbcp(i,j,k,isp) + phirmean(isp)
                if(discr_scheme.eq.'center')then
                   ur1phi(j,k)  = ubcp(i,j,k,1)*phirbcp(i,j,k,isp)
                endif
                ur2phi(j,k)  = ubcp(i,j,k,2)*phirbcp(i,j,k,isp)
                ur3phi(j,k)  = ubcp(i,j,k,3)*phirbcp(i,j,k,isp)
             end do
          end do
          ! Please note that phi in real space is stored in nlbc(0:nx+1,:,:,:)
          ! in this function. Later after finding the fluxes etc, in
          ! the same function, the real space phi is replaced into
          ! ubcp(0:nx+1,:,:,:). Please note the buffers into which I
          ! am sending the data in the two lines below this comment.
          ! Phirbc(2,:,:,:) is sent to ubcp(nx+2,:,:,:). These two
          ! extra buffers are needed for RPR operations. Since we
          ! will not be using ubcp(-1,j,k,:) and ubcp(nx+2,j,k,:) in this function,
          ! filling up those locations is safe. Also note a VERY
          ! important fact. phirbcp() is of size nx+2, whereas ubcp()
          ! is of size nx+4. So, while sending the non contiguous data blocks, these
          ! two arrays have different strides. Hence the use of
          ! VECSENDRECV2 which sends data between different
          ! stridetypes. VECSENDRECV sends data between similar stride
          ! types. See ../FLO/ibm.h for both these definitions.
          
          if(i.eq.2)then
             VECSENDRECV2(phirbcp(i,1,1,isp),1,twodrslice,fromproc,1,ubcp(nx+2,1,1,isp),1,urslice,toproc,1,decomp_group,status)
          else if(i.eq.nx-1)then 
             VECSENDRECV2(phirbcp(i,1,1,isp),1,twodrslice,toproc,0,ubcp(-1,1,1,isp),1,urslice,fromproc,0,decomp_group,status)
          end if
          call ff2rc(ur2phi, uf2phif)
          call ff2rc(ur3phi, uf3phif)
          do k  = 1, mz 
             do j = 1, my
                if((j.le.my2).and.(i.gt.0).and.(i.le.nx)) then 
                   nlphif(i,j,k,isp) = wy(j)*uf2phif(j,k) + wz(k)*uf3phif(j,k)
                end if
                if(discr_scheme.eq."center")then
                   fdderiv(i,j,k,isp) = ur1phi(j,k)
                else
                   fdderiv(i,j,k,isp)= zero
                end if
             end do
          end do
       end do !i=1, mx1
    end do !isp

#if !PARALLEL
    phirbcp(mx,:,:,1:nspmx) = phirbcp(1, :,:,1:nspmx)
#endif
    !if(discr_scheme.eq."center")
#if !PARALLEL
    if(discr_scheme.eq."center") fdderiv(mx,:,:,1:nspmx) = fdderiv(1&
         &,:,:,1:nspmx)
#endif
    if(discr_scheme.eq."upwind")then
!!$                         xi-1/2    xi+1/2          
!!$                     ......|..........|.......
!!$                     |          |            |
!!$                    xi-1        xi          xi+1   

!!$                   We need to compute (del f/del x)|i, where f = uT (T is the scalar)

!!$                      (del f/del x)|i = (f(i+1/2) - f(i-1/2))/dx
!!$ These are like finite volume cells centered around {xi}.
!!$ We cannot compute f(i+1/2) and f(i-1/2) directly since the polynomials are
!!$ discontinuous at these interfaces. 
!!$   For a second order accurate scheme, one must consider piecewise
!!$   linear polynomials.
#if PARALLEL
       if(.not.ALLOCATED(slope_u))then
!          BARRIER(decomp_group)
          ALLOCATE(uplus(0:nx), uminus(0:nx), U_plus(0:nx),&
               & U_minus(0:nx), slope_u(0:nx+1), slope_phi(0:nx+1))
       end if
#endif
       
       do isp = 1, nspmx
          do k = 1, mz
             do j = 1, my
#if PARALLEL
                
                slope_phi(0) = slope_factor*phi_ri*(phirbcp(1,j,k,isp)-ubcp(-1,j,k,isp))/(two*dx) 
                slope_phi(nx+1) = slope_factor*phi_ri*(ubcp(nx+2,j,k,isp)-phirbcp(nx,j,k,isp))/(two*dx) 
                IF(isp.eq.1)then
                   slope_u(0) = slope_factor*phi_ri*(ubcp(1,j,k,1)-uatminus1(j,k))/(two*dx)
                   slope_u(nx+1) = slope_factor*phi_ri*(uatnxp2(j,k)-ubcp(nx,j,k,1))/(two*dx)
                END IF
#endif
                do i = 1,nx !mx1
                   im1 = i-1
                   ip1 = i+1
#if !PARALLEL
                   if(im1.lt.1) im1 = mxf+im1-1
                   if(ip1.gt.mxf-1) ip1 = ip1-(mxf-1)
#endif
!!$ We have information at the grid locations {xi}. Consider piecewise linear polynomials in the intervals 
!!$   [xi-1/2,xi+1/2].

!!$ Form of the piecewise polynomial:
!!$   P[U](i) = U(i) + Sj*(x-x(i)); Sj is the slope of the  piecewise polynomial. Now how to choose the slope of the polynomial? 
!!$ First Order: Sj = 0

!!$ Second Order: 
                   slope_phi(i) = (phirbcp(ip1,j,k,isp)-phirbcp(im1,j,k,isp))/(two*dx) 
                   if(isp.eq.1)then
                      slope_u(i) = (ubcp(ip1,j,k,1)-ubcp(im1,j,k,1))/(two*dx) 
                   end if

                   if(limiter.eq."none")then

!!$                Sj = {T(i+1)-T(i-1)}/2dx --> second order slope 
                      phi_ri = one

                   else if(limiter.eq."minmod")then
                      
!!$              Sj = min[(T(i+1)-T(i))/dx, (T(i)-T(i-1))/dx] --> minmod limiter (for very high Re)
                      Dplus_vj = (phirbcp(ip1,j,k,isp)-phirbcp(i,j,k,isp))/dx
                      Dminus_vj = (phirbcp(i,j,k,isp)-phirbcp(im1,j,k,isp))/dx
                      slope_phi(i) = half*(sgn(Dplus_vj) + sgn(Dminus_vj))
                      slope_phi(i) = slope_phi(i) * MIN(ABS(Dplus_vj),ABS(Dminus_vj))
                      if(isp.eq.1)then
                         Dplus_vj = (ubcp(ip1,j,k,1)-ubcp(i,j,k,1))/dx
                         Dminus_vj = (ubcp(i,j,k,1)-ubcp(im1,j,k,1))/dx
                         slope_u(i) = half*(sgn(Dplus_vj) + sgn(Dminus_vj))
                         slope_u(i) = slope_u(i) * MIN(ABS(Dplus_vj),ABS(Dminus_vj))
                      end if

                   end if
                   slope_phi(i) = slope_factor*slope_phi(i)*phi_ri
                   if(isp.eq.1) slope_u(i) = slope_factor*slope_u(i)*phi_ri

                   if(isp.eq.1)then
!!$                From the cell i, one can get the following information at the interfaces (i-1/2) and (i+1/2)
!!$                (i-1/2) --> uplus stored at index i-1(im1)
!!$                (i+1/2) --> uminus stored at i
!!$                Uplus(xi+1/2) = P[U](i+1)|(xi+1/2) and
                      
!!$                Uminus(xi+1/2) = P[U](i)|(xi+1/2)
                      uplus(im1) = ubcp(i,j,k,1) - slope_u(i)*(dx/two) ! uplus(xi-1/2)
                      uminus(i) =  ubcp(i,j,k,1) + slope_u(i)*(dx/two) ! uminus(xi+1/2)
                   end if
                   U_plus(im1) = phirbcp(i,j,k,isp) - slope_phi(i)*(dx/two)  ! Uplus(i-1/2)
                   U_minus(i) =  phirbcp(i,j,k,isp) + slope_phi(i)*(dx/two) ! Uminus(i+1/2)
                   
!                   if(i.eq.1)then
!                      RSENDRECV(slope_u(i),1,fromproc,1,slope_u(nx+1),1,toproc,1,decomp_group,status) 
!                   else if(i.eq.nx)then 
!                      RSENDRECV(slope_u(i),1,toproc,0,slope_u(0),1,fromproc,0,decomp_group,status)       
!                   end if
                      
!!$ Evaluate the flux at an interface using the formula :
!!$           f(i+1/2) = 1/2[fplus(i+1/2) + fminus(i+1/2)] -1/2[ max(|uminus|,|uplus|)|(i+1/2){Uplus(i+1/2)-Uminus(i+1/2)}]

!!$           f(i-1/2) = 1/2[fplus(i-1/2) + fminus(i-1/2)] -1/2[ max(|uminus|,|uplus|)|(i-1/2){Uplus(i+1/2)-Uminus(i+1/2)}]

!!$ From the cell i, we can compute the following information:
!!$ At (i-1/2) --> fplus(i-1/2) stored at fdderiv(i-1) or fdderiv(im1)
!!$ At (i+1/2) --> fminus(i+1/2)stored at fdderiv(i) or fdderiv(im1)
                   
                   fdderiv(im1,j,k,isp) = fdderiv(im1,j,k,isp) + half*uplus(im1)*U_plus(im1)
                   fdderiv(i,j,k,isp) = fdderiv(i,j,k,isp) + half*uminus(i)*U_minus(i)

                   if(i.ne.1)then
!!$ Except for cell 1 (due to PBC), if I am at cell index i, then I am ensured that I now have the complete information 
!!$ for interface (i-1/2).
                      max_speed = MAX(ABS(uminus(im1)), ABS(uplus(im1)))
                      fdderiv(im1,j,k,isp) = fdderiv(im1,j,k,isp) - half*max_speed*(U_plus(im1)-U_minus(im1))
                   end if
                   
                end do
#if !PARALLEL
                max_speed = MAX(ABS(uminus(mx1)), ABS(uplus(mx1)))
                fdderiv(mx1,j,k,isp) = fdderiv(mx1,j,k,isp) - half*max_speed*(U_plus(mx1)-U_minus(mx1))
                fdderiv(mx,j,k,isp)  = fdderiv(1,j,k,isp)
#else
                uminus(0) = ubcp(0,j,k,1) + slope_u(0)*(dx/two)
                uplus(nx) = ubcp(nx+1,j,k,1) - slope_u(nx+1)*(dx/two)
                U_minus(0) = phirbcp(0,j,k,isp) + slope_phi(0)*(dx/two)
                U_plus(nx) = phirbcp(nx+1,j,k,isp) - slope_phi(nx+1)*(dx/two)
                fdderiv(0,j,k,isp) = fdderiv(0,j,k,isp) + half*uminus(0)*U_minus(0)
                max_speed = MAX(ABS(uminus(0)), ABS(uplus(0)))
                fdderiv(0,j,k,isp) = fdderiv(0,j,k,isp) - half*max_speed*(U_plus(0)-U_minus(0))
                fdderiv(nx,j,k,isp) = fdderiv(nx,j,k,isp) + half*uplus(nx)*U_plus(nx)
                max_speed = MAX(ABS(uminus(nx)), ABS(uplus(nx)))
                fdderiv(nx,j,k,isp) = fdderiv(nx,j,k,isp) - half*max_speed*(U_plus(nx)-U_minus(nx))
#endif

             end do
          end do
       end do
    end if

    DO isp = 1, nspmx
       DO i = 1, nx !mx1
          do k = 1, mz
             do j = 1, my 
                im1 = i-1
                ip1 = i+1
#if !PARALLEL
                if(im1.lt.1) im1 = mxf+im1-1
                if(ip1.gt.mxf-1) ip1 = ip1-(mxf-1)
#endif
                if(discr_scheme.eq."center")then
                   ur1phi(j,k) = (fdderiv(ip1,j,k,isp) - fdderiv(im1,j,k,isp))/(two*dx)
!!$                
                else if(discr_scheme.eq."upwind")then
                   ur1phi(j,k) = (fdderiv(i,j,k,isp) - fdderiv(im1,j,k,isp))/dx
                end if
                
             end do
          end do

          call ff2rc(ur1phi, uf1phif)

          nlphif(i,:,:,isp) = nlphif(i,:,:,isp) + uf1phif(:,:)
       end DO
#if !PARALLEL
       nlphif(mx, :,:,isp) = nlphif(1,:,:,isp)
#endif
    end DO
    phi_fluid_mean = zero 
    phi_solid_mean = zero 
    phimodmean = zero 
    phifmeanloc = zero
    phismeanloc = zero
    phimodmeanloc = zero
    DO isp = 1, nspmx
#if !PARALLEL
       DO i = 1, nx+1 !mx
#else
       DO i = 0, nx+1
#endif
          do k = 1, mz
             do j = 1, my

                ubcp(i,j,k,isp) = phirbcp(i,j,k,isp)
                if(i.gt.0.and.i.le.nx)then
                   if(j.le.my2) nlphif(i,j,k,isp) = -nlphif(i,j,k,isp)
                   if(fluid_atijk(i,j,k)) then 
                      phifmeanloc(isp) = phifmeanloc(isp) + ubcp(i,j,k,isp)
                   ELSE
                      phismeanloc(isp) = phismeanloc(isp)+ubcp(i,j,k,isp)
                   end if
                   phimodmeanloc(isp) = phimodmeanloc(isp) +  ubcp(i&
                        &,j,k,isp)
                end if
                !In nlphi, real transform of non linear term is stored in nlbc, therefore, here phirbc is stored in ubcp.
             end do
          end do
       end DO
#if !PARALLEL
       nlphif(mx,:,:,isp) = nlphif(1,:,:,isp)
#endif

       GLOBAL_DOUBLE_SUM(phifmeanloc(isp),phi_fluid_mean(isp),1,decomp_group)
       GLOBAL_DOUBLE_SUM(phismeanloc(isp),phi_solid_mean(isp),1,decomp_group)
    !   GLOBAL_DOUBLE_SUM(phimodmeanloc(isp),phimodmean(isp),1,decomp_group)
    end DO
    phi_fluid_mean(1:nspmx) = phi_fluid_mean(1:nspmx)/(count_fluid)
    phi_solid_mean(1:nspmx) = phi_solid_mean(1:nspmx)/(count_solid)
    phimodmean = phimodmean/(mx1*my*mz)
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A26,3(2x,g17.8))') "PHI MEAN F and SF = ",  phi_fluid_mean(1:nspmx), phi_fluid_mean(1:nspmx)*(one-maxvolfrac)
       WRITE(*,'(A26,3(2x,g17.8))') "PHI MEAN S and SF = ",  phi_solid_mean(1:nspmx), phi_solid_mean(1:nspmx)*maxvolfrac
       WRITE(*,'(A26,3(2x,g17.8))') "PHIMEAN OBTAINED", phi_fluid_mean*(one-maxvolfrac)+phi_solid_mean*maxvolfrac, phimean_des  
    end if

    RETURN 
    
  END SUBROUTINE form_nlphi


integer function sgn(a)
  implicit none
  REAL(prcn), INTENT(in) :: a
  if(a.lt.zero) then
     sgn = -1
  elseif(a.gt.zero)then
     sgn = 1
  else
     sgn = 0
  end if
end function sgn
   
END MODULE nl_allphi
