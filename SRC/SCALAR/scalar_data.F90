!GLOBAL DATA FOR SCALAR MODULES
!CURRENTLY SCALAR ROUTINES USE A LOT OF ARRAYS FROM THE HYDRODYNAMIC SIDE TO SAVE MEMORY.
!BE CAREFUL!
!AUTHOR: RAHUL GARG 
MODULE scalar_data
#include "../FLO/ibm.h"

  USE precision  
  USE constants 
  USE global_data,  ffphi=>ff
  IMPLICIT NONE 
  SAVE 
  
  INTEGER, PARAMETER :: nscalmax = 1 !maximum number of scalars possible
  REAL(PRCN), dimension(:), ALLOCATABLE :: sourcesink
  COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE ::  phif,nlphif,onlphif,sourcephif
  !COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE ::  ffphi
  COMPLEX(prcn), DIMENSION(:,:,:), ALLOCATABLE :: phiin,phiout
  !       Scalar real data 
  !REAL(prcn) :: phireal(mx, my, mz,nspmx)
  REAL(prcn), DIMENSION(:,:,:,:), ALLOCATABLE :: sourcephi
  !REAL(prcn), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gradphi
  !REAL(prcn) :: phi_an(mx,my,mz,nspmx)
  REAL(prcn), DIMENSION(:,:,:), ALLOCATABLE :: surf_scal_value, Nu3
  REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: phisurfall, phisurfall_nm1,phisurfall_nm2, &
          flux_body, flux_body2
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: phirmean,fphirmean,sum_flux_nm1,sum_flux_nm2,&
    sum_flux, GAMMA, phi_fluid_mean, phimodmean, phi_solid_mean
  
  REAL(PRCN), DIMENSION(:) , ALLOCATABLE :: nu_error_array, flux_global, flux_global2
  
  
  REAL(prcn) xis,xf,yi,yf,zi,zf,cons,uo,vo, nu_old, nu_error_hist, nu_error, phimean_des

  INTEGER ::  nxi,nxf,nyi,nyf,nzi,nzf,diff,nspmx, iter_scal
  REAL(prcn) :: phisurf,phistream,  Pr_or_Sc_nu!,gamma(nspmx)
  LOGICAL :: sourcepresent, LUMPED_CAP, zero_flow_in_solid, dorpr_scal, setphimean
  
  
  !Integer :: inn_itn_step

  NAMELIST/scal_propt/ nspmx,phistream, setphimean, phisurf,sourcepresent, LUMPED_CAP, &
    Pr_or_Sc_nu, zero_flow_in_solid, dorpr_scal

END MODULE scalar_data

