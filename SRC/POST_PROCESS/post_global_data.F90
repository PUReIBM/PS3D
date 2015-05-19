MODULE post_global_data
#include "../FLO/ibm.h"
  USE precision 
  USE constants 
  USe global_data
  USE scalar_data 
  IMPLICIT NONE
  SAVE 
  PUBLIC

  INTEGER :: nmis,nmismax,percmax, mxmax, ndim1,ndim2,nphsmax
  PARAMETER(nmismax=30,percmax=11, mxmax = 210,ndim1=ndim-1,nphsmax = 10)
  parameter(ndim2=ndim*ndim)
  CHARACTER*8 :: POST_RUNNAME, PER_CONF
  ! Freely Evolving Suspensions
  REAL(prcn) :: movestart_time, post_proc_time
  INTEGER, SAVE :: grantempunit=1, aivjunit=1, anisunit=1,&
       & source_dissip_unit=1, eigenvalue_unit=1, gofrunit=1,&
       & source_dissip_pres_visc_unit=1, accphiunit=1

  Type :: REALIZATION
     INTEGER :: mis, nphases
     TYPE(PHASE), DIMENSION(:),POINTER :: phase_info
     REAL(prcn),DIMENSION(:,:),POINTER :: for,vel
     REAL(prcn) :: usmean(ndim), ufmean(ndim), mixmeanslip(ndim), ufvar(ndim), uiuj(ndim,ndim)
     REAL(prcn) :: gran_temp, maxvolfrac
     REAL(prcn) :: mean_force(ndim), sigma_F(ndim) ! Standard deviation of each particles drag
     REAL(prcn) :: mean_pres(ndim), sigma_pres(ndim) ! Standard deviation of each particles drag
     REAL(prcn) :: mean_visc(ndim), sigma_visc(ndim) ! Standard deviation of each particles drag
     TYPE(REALIZATION), POINTER :: next
  END Type REALIZATION

  Type :: ENSEMBLE
     INTEGER :: nmis
     INTEGER :: nsamples
     REAL(prcn) :: mean_force(ndim),sigma_F(ndim) ! Standard deviation in the mean drag
     REAL(prcn) :: mean_pres(ndim), sigma_pres(ndim) ! Standard deviation of each particles drag
     REAL(prcn) :: mean_visc(ndim), sigma_visc(ndim) ! Standard deviation of each particles drag
     REAL(prcn) :: meanslip(ndim), meanslipmod,grantemp, volfrac

     REAL(prcn) :: ufvar(ndim), sigma_ufvar(ndim)
     REAL(prcn) :: norm_drag
     REAL(prcn) :: force_max(ndim), force_min(ndim), vel_max(ndim), vel_min(ndim)
  END Type ENSEMBLE


  TYPE(ENSEMBLE) :: ENSAVG
  TYPE(REALIZATION), POINTER :: mis_data
  TYPE(REALIZATION),POINTER :: current_mis

  REAL(prcn) :: tke_fluid, tke_fluid_var

  REAL(prcn) :: nusselt_all_mis(nmismax), nusselt_perp(nmismax), nusselt_para(nmismax), nlcontrib(nmismax, nscalmax), diffcontrib(nmismax, nscalmax), fphicontrib(nmismax, nscalmax),  fphicontrib_solid(nmismax, nscalmax)
  REAL(prcn) :: scal_frac_nl(mxmax, 3, nmismax), scal_frac_diff(mxmax, 3, nmismax), solid_area_frac(mxmax, nmismax), scal_frac_nl_diff(mxmax,nscalmax,nmismax), scal_frac_fphi(mxmax,nscalmax,nmismax)

  REAL(prcn) :: mis_norm_drag_spec(nmismax,nphsmax),mis_norm_drag_chem_spec(nmismax,nphsmax)
  REAL(prcn) :: avg_norm_drag_spec(nphsmax),avg_norm_drag_chem_spec(nphsmax)
  REAL(prcn) :: norm_drag_spec_errbar(nphsmax),norm_drag_chem_spec_errbar(nphsmax)

  INTEGER :: MIS_CONV(nmismax), convmiscount
  Type :: conf_interval
     character*8 :: confper(percmax)
     real(prcn) :: confi(percmax)
  END Type conf_interval

  TYPE(conf_interval) :: cis(nmismax)

  REAL(prcn) :: confin

  INTEGER :: nmeasvols
  PARAMETER(nmeasvols=5)

  REAL(prcn) ::  NUMBER_IN_MIS(nmismax,nmeasvols),&
       & NUMBERSQ_IN_MIS(nmismax,nmeasvols),&
       & VOLUME_IN_MIS(nmismax,nmeasvols),&
       & VOLUMESQ_IN_MIS(nmismax,nmeasvols)

  REAL(prcn) ::  NUMBER_IN_MEASVOL(nmeasvols),&
       & NUMBERSQ_IN_MEASVOL(nmeasvols),&
       & VOLUME_IN_MEASVOL(nmeasvols),&
       & VOLUMESQ_IN_MEASVOL(nmeasvols)
  REAL(prcn) :: NUMBER_STDDEV(nmeasvols), NUMBERSQ_STDDEV(nmeasvols),&
       & VOLUMESQ_STDDEV(nmeasvols), VOLUME_STDDEV(nmeasvols), LMEAS_VOLS(nmeasvols)

  REAL(prcn) :: VEL_PAR_MEASVOL(nmeasvols),&
       & VEL_PER_MEASVOL(nmeasvols), VEL_PAR_STDDEV(nmeasvols),&
       & VEL_PER_STDDEV(nmeasvols)

END MODULE post_global_data
