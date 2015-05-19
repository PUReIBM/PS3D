MODULE nlphi
#include "../FLO/ibm.h"
  Use precision
  USE scalar_data
  USE fftw_interface
!  USE global_data! ONLY : mx, my, mz, my2, mxf, foffset, aliasflag, nx
  USE nl_allphi
  USE bc_scalar
  !USE bc_scalar_new
  
  Use nlmainarrays, Only : onlphirbcp=>onlbc, nlphirbcp=>nlbc, phirbcp=>ubc
  !Use nlphimainarrays, Only : onlphip,nlphip, phip
  IMPLICIT NONE 
  !--------------------------------------------------------------------

  !calculate the nonlinear (convective) terms in the Scalar equation
  !Rahul Garg 
  !	
  !use phase-shifting to dealias
  !-----------------------------------------------------------------------
CONTAINS 
  SUBROUTINE nlphi_graduphi(rks)


    IMPLICIT NONE 
    INTEGER, INTENT(in) ::  rks

    REAL(prcn) ::  nlphimean(nspmx),nlphimeanloc(nspmx)
    REAL(prcn) ::  eold,fr1(my,mz)
    COMPLEX(prcn) ::  tmpc!, onlf(my2,mz)
    INTEGER i,j,k,l,ii,n,isp


    DO isp=1,nspmx
       DO k=1,mz
          DO j=1,my2
             DO i=1,nx+1 !mx
#if PARALLEL
                if(i.le.nx)then
#endif
                   onlphif(i,j,k,isp)=nlphif(i,j,k,isp)
                   nlphif(i,j,k,isp)=czero
#if PARALLEL
                end if
#endif
             ENDDO
          ENDDO
       ENDDO

    ENDDO
    IF(Re.EQ.ZERO) THEN 
       WRITE(*,'(A50)') 'TRANSFORM ONLY PHI ARRAY: NOT CALCULATING NLPHI'
       do isp = 1, nspmx
#if !PARALLEL
          do i = 1, nx !mx1
#else
          do i = 0,nx+1
#endif
             call ff2cr(phif(i,:,:,isp), phirbcp(i,:,:,isp))
             do k = 1, mz
                do j = 1, my
                   nlphirbcp(i,j,k,isp) = zero
                   onlphirbcp(i,j,k,isp) = zero
                end do
             end do
             if(i.eq.2)then
                VECSENDRECV(phirbcp(i,1,1,isp),1,urslice,fromproc,1,phirbcp(nx+2,1,1,isp),1,toproc,1,decomp_group,status)
             else if(i.eq.nx-1)then 
                VECSENDRECV(phirbcp(i,1,1,isp),1,urslice,toproc,0,phirbcp(-1,1,1,isp),1,fromproc,0,decomp_group,status)
             end if

          end do
          
       end do

    ELSE

       CALL form_nlphi
       DO isp=1,nspmx
#if PARALLEL       
          nlphirbcp(0,:,:,isp) = zero
          onlphirbcp(0,:,:,isp) = zero
#endif
          
          DO i=1,nx !mxf
             CALL ff2cr(nlphif(i,:,:,isp), nlphirbcp(i,:,:,isp))
             CALL ff2cr(onlphif(i,:,:,isp),onlphirbcp(i,:,:,isp))
          ENDDO
#if PARALLEL       
          nlphirbcp(nx+1,:,:,isp) = zero
          onlphirbcp(nx+1,:,:,isp) = zero
#else
          nlphirbcp(nx+1,:,:,isp) = nlphirbcp(1,:,:,isp)
          onlphirbcp(nx+1,:,:,isp) = onlphirbcp(1,:,:,isp)
#endif
       ENDDO

       do isp = 1, nspmx
          nlphimeanloc(isp) = SUM(nlphirbcp(1:nx,:,:,isp))
          GLOBAL_DOUBLE_SUM(nlphimeanloc(isp),nlphimean(isp),1,decomp_group)
       end do
       nlphimean(:) = nlphimean(:)/(mx1*my*mz) 
       if(I_AM_NODE_ZERO)Write(*,'(A25,10(2x,g12.5))')'NLPHIMEAN = ',( nlphimean(isp), isp = 1,nspmx)
       
    end IF


    CALL bcsetscal(rks)

  END SUBROUTINE nlphi_graduphi
END MODULE nlphi
    



