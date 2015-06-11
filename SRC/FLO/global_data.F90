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


MODULE global_data
#include "ibm.h"
	USE precision 
	USE constants 
	!USe reversal_pts
	IMPLICIT NONE
#if PARALLEL
	include "mpif.h"
#endif
!#if FFTW3
!  include "../../include/fftw3.h"
!#endif
	SAVE 
	PUBLIC



	INTEGER ::  mx,my,mz !,mx1,mx2,my2,mxf, nx
	INTEGER ::  ndim
	INTEGER ::  n1,n2
	!INTEGER ::  xstart, xend
	PARAMETER(ndim=3)

	PARAMETER(n1=45000,n2=20000)

	!-----------------------------------------------------!
	!All complex Declarations 
	!-----------------------------------------------------!
	!DOUBLE COMPLEX TYPE DECLARATIONS

	COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: u, ff
	COMPLEX(prcn), DIMENSION(:,:,:), ALLOCATABLE :: p, pstar !uin, uout
	COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE , TARGET :: nl,onl
	COMPLEX(prcn), DIMENSION(:), ALLOCATABLE  :: wx,wy,wz, shifty, shiftz !, wx2, wy2, wz2
	REAL(prcn), DIMENSION(:,:,:), ALLOCATABLE ::  w2
	COMPLEX(prcn), DIMENSION(:,:), ALLOCATABLE :: shiftyz
	COMPLEX(prcn) :: total_mean_forcing(ndim), total_div,divloc
	!complex(prcn) :: utin(mx,my,mz,ndim)
	!complex(prcn) :: divglobal(mx,my2,mz,ndim)
	!complex(prcn) :: phif(mx,my2,mz,1)
	LOGICAL, Dimension(:,:,:), ALLOCATABLE :: fluid_atijk
	INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE   :: GNACC_PART
	INTEGER :: MAXGN_PARTNEIGH
	PARAMETER(MAXGN_PARTNEIGH = 8)

!!$  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE   :: wtbar
!!$  INTEGER, DIMENSION(:,:,:), ALLOCATABLE   ::  wtbar_tmp

	!Timing Related 
	INTEGER :: ID_HOUR, ID_MINUTE, ID_SECOND, ID_YEAR, ID_MONTH, ID_DAY
	CHARACTER*64 :: ID_NODE



	!REAL(prcn), DIMENSION(:), ALLOCATABLE :: g
	REAL(prcn) :: dt_tmp_conv, dt_tmp_diff, dt_tmp_vis,dt_tmp_grav, dtcoll
	REAL(prcn) :: maxvolfrac, upi_sf, ferror_hist, nl_slope_factor
	!Time steps 
	REAL(prcn) :: dtc, dtv, dtd!,  rnp(mxf,nprocy, nproczmy,mz,2) 
	!Time scales

	real(prcn) :: dt_old, AB_relax !relaxation factor for Adam Bashforth 0 means no modification

	REAL(prcn) :: t_conv, t_vis, t_diff, t_min, t_turb, dtnew, dt_tmp, t_grav, AUTO_CORR_SEP_TIME, rhoofs0(ndim)

	REAL(prcn) :: dem_time, dem_time_current, hydro_time, hydro_time_current, &
				&	nl_time, nl_time_current, bc_time, bc_time_current, vp_time, &
				&	vp_time_current, scalar_time, scalar_time_current, hydroforce_time, &
				&	hydroforce_time_current, fourier_time, fourier_time_current, fourier_time_nl, &
				&	new_timestep_time, new_timestep_time_current, mpg_time, pgrad_time, reversal_time

	REAL(prcn) ::  coef(3,4) !,coeff((my+15)+2*(mz+15))
	REAL(prcn) ::  umean(ndim), umax, dtorig, tsprgmin, umodmean,&
				& frmean(ndim), frmeanfluid(ndim),pres_total(ndim),&
				& visc_total(ndim), visc_total_old(ndim),  pres_avg(ndim),&
				& visc_avg(ndim), pres_drag, visc_drag, total_drag
	REAL(prcn) ::  ax(ndim), cf , icf, cforig
	REAL(prcn) ::  dx2, bndrad, doml(3)
	REAL(prcn) :: beta,wmax2, mpg_avg, total_drag_mpg
	!real(prcn) ::  norm(ndim,nsmx),mask(mxf,my,mz)
	!real(prcn) ::  f(nbmx,ndim,nsmx),intf(nbmx,ndim,nsmx)
	!real(prcn) ::  fint(nbmx,ndim,nsmx),ifint(nbmx,ndim,nsmx)
	REAL(prcn) :: urf,prf,  resid_flow!, mean_pressure(mx1)

	LOGICAL :: FIRST_PASS, VEL_CONVERTED, DEBUG_CHECK, killjob, check_point
	REAL(prcn) ::  mtot, wp(ndim), apmax,ucm(ndim), rpart, uph(ndim),  x
	REAL(prcn) ::  tau(2*ndim), energy, upbcset

	INTEGER(8) :: count_solid, count_fluid, count_total  

	REAL(prcn), dimension(:,:), Allocatable ::  force,pres,visc, cd,&
				& torq,force_loc, presloc,viscloc,torqloc, force_chem, contact_force 
	REAL(prcn), Allocatable ::  visc_force_bnd(:,:,:), pres_force_bnd(:,:,:), xc_bnd(:,:,:)

	REAL(prcn), dimension(:,:,:), Allocatable :: ap, up, angv,anga

	REAL(prcn), dimension(:), Allocatable ::  cs,mp, mpart,mompart
	REAL(prcn), dimension(:,:), Allocatable  ::  xs, xc, xc_init, tmpfor, velbdy, velbdy0, color(:)
	REAL(prcn), DIMENSION(:), ALLOCATABLE :: radbdy, radibdy, radobdy, rado2bdy

	REAL(prcn), ALLOCATABLE, DIMENSION(:,:) :: XC_GENER
	REAL(prcn), ALLOCATABLE, DIMENSION(:) :: RAD_GENER
	REAL(prcn) :: tphase,  ramp_frac_time

	LOGICAL ::  gauss_u, gauss_p, gauss_phi, diff_ts, gauss_scal, frombnd, fromrpr

	LOGICAL ::  movingcv


	logical, allocatable :: myid_particles(:)

	!ALL INTEGER DECLARATIONS 

	INTEGER:: ounit, eunit, runit, ramp_frac_steps, unit_overlap
	INTEGER ::  nbnd,nrpr
	INTEGER ::  cuty,cutz
	!integer ::  nbin(my2,mz),nmode(my2)

	!integer ::  nn(my,my2,my)
	INTEGER ::  lc, agflag, idumstep, iglobstep, itrmax, iter_gauss !, iter_u(3), iter_p, iter_phi

	INTEGER, DIMENSION(:), ALLOCATABLE:: revp
	INTEGER:: revgroup

	!INTEGER*8 :: planr2c_2d, planc2r_2d, planr2c_3d, planc2r_3d

	!Type(parrayrev), Dimension(:), Allocatable:: pgrev


	real(prcn), allocatable :: gofr_avg(:),gofr_mis(:,:), rad_bin(:), int_dist(:)
	real(prcn) :: matern_treshold
	integer :: nbins, nsmall_length


	!-------
	! Namelist input variables
	!-------
	CHARACTER(LEN=40):: RUN_NAME, errfile,outputfile
	INTEGER :: minunitno, maxunitno, saveitns, mbuffer
	!Flags
	INTEGER :: aliasflag, maxzeroflag, revernorm, iplot,imove, &
		& iscalon, irestart, iscal_restart, icomplex, &
		& saveforrestart
	LOGICAL ::  xperiodic(ndim-1), adaptive, ibidisperse,set_umean, set_mpg,rk3, &
		& mean_vel_to_particles, move_particles, dorpr, dobnd,&
		& gof_avg, write_snapshot, only_dem, iturbon, include_nmis, pdf_normalized, bubble_particles, &
		& init_parts_with_fluid_vel
	CHARACTER(LEN=11):: rstsave, outformat, rstread
	!Part_propt
	REAL(prcn) :: mparti, momparti, radm, gap1, r,dr,f1,f2, tstartmove, &
		& saveforres_time,sprgconst, lybyd, dbydx, dia_phys, dchar,&
		& Tratio, voldom , hbydx, hbyd, widthbyd, aspect_ratio, volume_fraction
	INTEGER :: nbody, totcolls, count_restart, inneritermax
	CHARACTER(LEN=80) :: input_type, simtype, discr_scheme, limiter, collision_type, vel_distrib

	!thermal quantities for the bodies(sub s)  and fluid (sub f)
	REAL(prcn) rhos,cps,ks, rhof, cpf, kf
	!floandgrid

	REAL(prcn) :: dt, t, tstart, tend, tendused, TOL_FERROR,&
		& TOL_GRAN_ERROR, wtime_maxhrs, vis, Re, re_out, dx, dy,dz, upi(ndim),&
		& cfl, vfl, ag(ndim), mpg(ndim), uchar(ndim),&
		& gran_temp, grant_i, ReT, uconv, ucharmod, yalpha(2), GRAV(ndim)

	REAL(prcn) :: ufmean(ndim),usmean(ndim), meanslip(ndim),meanslipmod&
		&, usmean_des(ndim), usmean_act(ndim), ufmean_des(ndim),&
		& ufmean_old(ndim), fsslip(ndim), fsslipmod, mesh_vel(ndim), umeanslip, slip_des(ndim)
	REAL(prcn) :: frame_accln(ndim), frame_vel(ndim), frame_pos(ndim)
	REAL(prcn) :: mixmeanslipmod,dufmeandt(ndim)
	LOGICAL :: impose_grav, include_frmeanfluid, cage_simul,&
		& fixed_moving, moving_from_fixed_bed, use_drag_law, communicator_done, use_fes, zero_slip
	REAL(prcn) :: mixmean_vel(ndim)
!  INTEGER :: foffset, foffsetorig
  
	Integer, SAVE :: unitdrag=1, unitdragsum=1, unitnormdrag=1,&
		& unitnormdragchem=1, unitdragtavg=1, unitdragz=1, unitforce=1&
		&, unitdrag_comps=1, unitnormdragpoly = 1, unitgrantemp=1,&
		& unitpartinfo=1, sphrunit = 1, partunit = 1

	REAL(prcn) :: flo_ang(3), archno
	!Outptctrl
	INTEGER ::  lcinput,stepbody
	REAL(prcn) :: tstartplot, ferror, fold
	!Interpolation related quantities
	INTEGER:: order, ob2l, ob2r
	REAL(prcn), DIMENSION(:,:,:,:), ALLOCATABLE:: gstencil, vstencil
	CHARACTER(LEN=7):: scheme, coupling, interp_scheme

	REAL(prcn), DIMENSION(:,:,:), POINTER :: weightp
	REAL(prcn), DIMENSION(:,:,:,:), ALLOCATABLE :: vsten, nlsten&
		&,onlsten, dfsten, ppgrsten , sourcesten
	REAL(prcn), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gradphisten
	REAL(prcn), DIMENSION(:,:,:), ALLOCATABLE ::  prsten
	REAL(PRCN), DIMENSION(:) , ALLOCATABLE :: ferror_array
	lOGICAL :: inner_itn, final_restart, nlcons
	LOGICAL:: intx_per, inty_per, intz_per, write_output,&
		& double_delta, heaviside, scalneumann, change_gamma, dtmin,&
		& heaviside3, smear_normal, ibordel, equal_momentum
	INTEGER :: lc_steps, mis_hs, mis_mat, nerr_steps
	REAL(prcn) :: delta, spread, clip, flux_des, hhat, min_part_sep

	REAL(prcn) :: force_factor


	!event driven algo related 
	LOGICAL :: flow_converged,scal_converged, soft_start&
		&,base_state_saved,compute_auto_corr

	LOGICAL :: FROM_POST, POST_NO_FLOW_MEM_ALLOC, CHANGE_FRAME

	REAL(prcn) :: coeff_rest, pvel_var, volfrac_rat,&
		& dia_ratio, vsim, RED1
	INTEGER :: min_colisions
	REAL(prcn) :: source_hydro(ndim)
	! PSD related 
	CHARACTER*80 :: psd_type
	REAL(prcn) :: sigmabydavg,phiavg
	INTEGER :: nphases, nphsc2
	REAL(prcn) :: char_length
	TYPE PHASE
		REAL(prcn) :: dia, volfrac, volfracg
		REAL(prcn), DIMENSION(ndim) :: mean_spec_vel
		INTEGER :: npart, pstart, pend, nbnd, nrpr
		REAL(prcn) :: ferror, fold, ferror_hist
		REAL(prcn) :: grant_error, grant_old, gran_error_hist
		REAL(PRCN), DIMENSION(:) , POINTER :: ferror_array
		REAL(PRCN), DIMENSION(:) , POINTER :: grant_array
		REAL(prcn), DIMENSION(:,:), POINTER :: bndpts
	END TYPE PHASE
	TYPE(PHASE), DIMENSION(:), ALLOCATABLE, TARGET :: phase_array
  
	TYPE PARTICLE
		INTEGER :: iphs, nrpr_active, drag_active
		LOGICAL(KIND=1),DIMENSION(:), POINTER :: if_rev, if_drag
	END TYPE PARTICLE
	TYPE(PARTICLE), DIMENSION(:), ALLOCATABLE, TARGET :: part_array

	REAL(prcn), DIMENSION(:), ALLOCATABLE:: percent_buf
	REAL(prcn), DIMENSION(:), ALLOCATABLE:: norm_drag_spec,norm_drag_chem_spec
	REAL(prcn), DIMENSION(:,:), pointer :: bndarray
	REAL(prcn) :: mean_volfrac, norm_drag, norm_drag_chem
	!-------
	! PARALLEL Declarations
	!-------
	INTEGER :: myid, nproc, nprocy, nprocz
	integer, allocatable :: nproc_array(:)

#if PARALLEL
	INTEGER :: err_code, comm_group, decomp_group, comm_cart_2d
	parameter(comm_group = MPI_COMM_WORLD)
	INTEGER :: node_zero, toy, toz, fromy, fromz
	parameter(node_zero=0)
	INTEGER :: xy_logical, xz_logical, xy_double, xZ_double, xy_vel, xz_vel
	INTEGER :: status(MPI_STATUS_SIZE)
#endif
	integer :: nmis, skip_num, sa_max_itr, igeometry, isa, sa_in_itr, icohesive
	real(prcn) :: rad_factor, power, rad_reduction_ratio

	INTEGER :: count_phi, count_theta, igoemetry, cluster_part_min, cluster_part_max, slice_dx
	REAL(prcn), Dimension(:), Allocatable :: phi_array, theta_array
	REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: Fvis_theta, Fvis_phi, pres_theta, pres_phi 

	real(prcn) :: post_error, post_small, rad_ratio, func_rad_ratio, cooling_rate, &
			& temp_init, sa_error, objective_coef1, objective_coef2, lybyd_small, d0_rat, hamaker

	logical :: rotated
	real(prcn) :: tke, epsf, tke_i, epsf_i

  !-------
  ! Namelist definitions
  !-------
	NAMELIST/Flags/ aliasflag, maxzeroflag, revernorm, imove&
		&,iscalon, irestart, iscal_restart, saveforrestart&
		&,saveforres_time, adaptive, ab_relax, set_umean, set_mpg,nlcons, rk3&
		&, mean_vel_to_particles, dorpr, dobnd, include_frmeanfluid,&
		& impose_grav, debug_check, cage_simul, fixed_moving, gof_avg,&
		& only_dem, movingcv, iturbon, moving_from_fixed_bed, from_post, &
		& change_frame, use_drag_law, use_fes, simtype, zero_slip, init_parts_with_fluid_vel

	NAMELIST/Files/outputfile, rstsave, outformat, saveitns, minunitno, &
		& maxunitno, errfile, rstread, write_output, write_snapshot


	NAMELIST/Part_propt/ nbody,dr, f1, f2,&
		& tstartmove, lybyd, dbydx, dia_phys, rhos, cps, ks,&
		& input_type, discr_scheme, limiter, min_part_sep&
		&,nphases,collision_type, vel_distrib, nl_slope_factor, bubble_particles

	NAMELIST/floandgrid/ mx, my, mz, nprocy, nprocz, tstart,ramp_frac_steps, tend& !, mxf
		&,nerr_steps, wtime_maxhrs, saveforres_time, TOL_FERROR,  vis,&
		& Re, rhof, cpf, kf, archno, &
		& cfl, vfl, mbuffer, mpg, ReT, flo_ang,&
		& TOL_GRAN_ERROR
  
!!$  Namelist/testing/ oldmethod, onlydiff, creepflow, heaviside3, double_delta, inner_itn,&
!!$       & heaviside, scalneumann, change_gamma, smear_normal, dtmin,&
!!$       & delta, newmethod, spread, clip, flux_des
	NAMELIST/psd/  psd_type,nphases,sigmabydavg,phiavg,equal_momentum
	NAMELIST/gcg/  coeff_rest, pvel_var, min_colisions, mis_hs
	NAMELIST/bidisperse/  volfrac_rat, RED1, dia_ratio, yalpha
	NAMELIST/lubtest/  hbyd
	NAMELIST/matinput/ hhat, mis_mat, ibordel, matern_treshold
	NAMELIST/riser/  volume_fraction, widthbyd, aspect_ratio

	namelist/postlist/ nmis, nbins, skip_num, post_no_flow_mem_alloc, &
		&	rad_factor, power, sa_max_itr, post_error, &
		&	post_small, rad_ratio, func_rad_ratio, cooling_rate, temp_init, sa_error, objective_coef1, &
		&	objective_coef2, igeometry, isa, sa_in_itr, rad_reduction_ratio, icohesive, nsmall_length, &
		&	lybyd_small, cluster_part_min, cluster_part_max, rotated, include_nmis, pdf_normalized
END MODULE global_data









