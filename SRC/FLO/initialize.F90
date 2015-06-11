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


module initialize
#include "ibm.h"
	USE precision  
	Use constants 
	use global_data
	Use errormesgs
	Use general_funcs
	!USE scalar_data
	use init_turb
	implicit none 
contains
	subroutine init_params
		implicit none 
		call set_default

#if PARALLEL
		status= MPI_STATUS_IGNORE
#endif
		call read_nmls 

		!iscal_restart = irestart 
		!if (I_AM_NODE_ZERO) then
		!   call write_nmls
		!end if
		!--------
		! write namelists into output file
		!--------
		!CALL write_input(ounit)
	end SUBROUTINE init_params
  
	SUBROUTINE set_default 
		Implicit None
		INTEGER :: isp
		!Flags
		mx = UNDEFINED_I
		my = UNDEFINED_I
		mz = UNDEFINED_I
		nphases = 1
		ramp_frac_steps = 20.d0
		nerr_steps = 100
		wtime_maxhrs = 23.d0 

		TOL_FERROR =    1.E-06 
		TOL_GRAN_ERROR = 1.E-03
		aliasflag = 0
		maxzeroflag = 0
		revernorm = 1  
		dorpr = .TRUE.
		dobnd = .TRUE.
		imove = 0
		iscalon = 0
		irestart = 0
		saveforrestart       = 1
		saveforres_time = 6.d0
		xperiodic = .TRUE.
		adaptive = .TRUE.
		nlcons = .TRUE.
		set_umean = .FALSE.
		set_mpg = .FALSE.
		rk3 = .FALSE.
		mean_vel_to_particles = .FALSE.
		include_frmeanfluid = .FALSE.
		impose_grav = .FALSE.
		cage_simul = .FALSE.
		gof_avg = .FALSE.
		write_snapshot = .FALSE.
		fixed_moving = .FALSE.
		only_dem = .FALSE.
		movingcv = .FALSE.
		debug_check = .TRUE.

		!Files
		outputfile = "ffd.out"
		rstsave ="unformatted"
		outformat     = "formatted"
		errfile       = "ffd.err"
		rstread       = "unformatted"

		saveitns      = 15000
		write_output = .false.
		!outfilesize   = 250           ! In Megabytes
		minunitno     = 120
		maxunitno     = 300

		!Outputctrl/

		!Part_propt/
		nbody = 1 
		dr = 1.0
		f1 = 1.6d0
		f2 = 1.6d0
		dia_phys = 1.d0
		tstartmove = 0.00744
		lybyd = UNDEFINED_R
		dbydx = UNDEFINED_R
		rhos = 1.0
		cps = 447.0
		archno = 20.0
		ks = 0.03
		phiavg = 0.1d0
		volfrac_rat = 1.0
		nl_slope_factor = 1.d0 
		input_type = "default"
		discr_scheme = "center"
		limiter = "none"
		vel_distrib = "maxwellian"
		collision_type = "none"
		psd_type = "discrete"
		widthbyd = lybyd
		aspect_ratio = one
		volume_fraction = 0.05
		min_part_sep = 0.d0
		dia_ratio = 1.d0

		!floandgrid/ 
		vis = 0.012
		Re = 100.00
		cfl= 0.5
		vfl = 1.0
		mbuffer = 0
		ReT = zero
		mpg = 0.0

		flo_ang(1) = 0.d0

		flo_ang(2) = 90.d0
		flo_ang(3) = 90.d0
!!$    !TESTING
!!$    OLDMETHOD = .False.
!!$    onlydiff = .False.
!!$    heaviside3 = .True.	
!!$    double_delta  = .True.
!!$    inner_itn  = .True.
!!$    heaviside = .False.
!!$    scalneumann = .False.
!!$    change_gamma = .False.
!!$    dtmin = .True.
!!$    delta  = 0.12

		!gcg
		ibidisperse  = .FALSE.
		coeff_rest = 1.0
		pvel_var = 1.0
		min_colisions = 10
		mis_hs = 1 
		!mat
		hhat = 2.d0
		mis_mat = 1
		ibordel = .false.
		matern_treshold = 0.9

#if 0
		!scal_part_propt
		nspmx=1
		phistream = 1.d0
		phisurf = 0.d0
		sourcepresent = .FALSE.

		PR_OR_SC_NU =    0.700000000000000 
		zero_flow_in_solid = .FALSE.
		dorpr_scal = .false.
		setphimean = .true.
#endif

		!^^^^^^^ 06-24-2010 Mohammad
		post_no_flow_mem_alloc = .false.
		from_post = .false.	
		iturbon = .false.
		!tke_converged = .false.
		moving_from_fixed_bed = .false.
		nmis = 1
		rad_factor = 2.0
		power = 1.0
		nbins = 200
		sa_max_itr = 1
		post_error = 1d-3
		post_small = 1d-3
		rad_ratio = one
		cooling_rate = 0.95

		func_rad_ratio = 6
		temp_init = 1000
		objective_coef1 = one
		objective_coef2 = one
		igeometry = 0
		isa = 0
		sa_in_itr = 1
		rad_reduction_ratio = 0.9
		icohesive = 0
		nsmall_length = 1
		cluster_part_min = 0
		cluster_part_max = 0
		lybyd_small = one
		rotated = .false.
		hamaker = 1d-19
		d0_rat = 1d-4

		slice_dx = 1
		use_drag_law = .false.
		force_factor = one
		use_fes = .false.
		inneritermax = 10
		simtype = 'accframe'
		zero_slip = .false.
		nproc = 1
		nprocy = 1
		nprocz = 1
		bubble_particles = .false.
		init_parts_with_fluid_vel = .false.


		epsf_forcing = -10.0
		include_nmis = .false.
		grant_i = zero

		ab_relax = zero
		pdf_normalized = .true.
		skip_num=1
		!------------------------------------------------------------------------

		dem_time = zero
		dem_time_current = zero
		hydro_time = zero
		hydro_time_current = zero
		nl_time = zero
		nl_time_current = zero
		bc_time = zero
		bc_time_current = zero
		vp_time = zero
		vp_time_current = zero
		new_timestep_time = zero
		new_timestep_time_current = zero
		!scalar_time = zero
		!scalar_time_current = zero
		hydroforce_time = zero
		hydroforce_time_current = zero
		fourier_time = zero
		fourier_time_current = zero
		fourier_time_nl = zero
		!------------------------------------------------------------------------
	end SUBROUTINE set_default


	SUBROUTINE read_nmls
		implicit none 
		Logical:: nmlexist
		Integer:: unitno, ierr
		CHARACTER*80 :: FILENAME

		unitno = 101     ! have to use arbitrary unit number

		if (TRIM(RUN_NAME).EQ."")then
			CALL GENER_FILENAME(FILENAME,'floparam.in')
		else
			CALL GENER_FILENAME(FILENAME,TRIM(RUN_NAME)//'_floparam.in')
		end if
		open(unit=unitno,file=FILENAME,form="formatted",DELIM='APOSTROPHE', IOSTAT=ierr)
		read(unitno,NML=Floandgrid)
		read(unitno,NML=Part_propt)
		read(unitno,NML=Flags)
		if(imove.ne.1)cage_simul = .FALSE.
		if(imove.ne.1)fixed_moving = .FALSE.
		if(set_umean) set_mpg = set_umean
		read(unitno,NML=Files)
		read(unitno,NML=psd)     
		read(unitno,NML=gcg)       
		if ((TRIM(input_type).eq.'lubtest')) then
			read (unitno,NML=lubtest)
			read (unitno,NML=bidisperse)
		endif
		if((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'bidisp'))then
			read (unitno,NML=bidisperse)
		end if
		if(TRIM(input_type).eq.'mat'.or.input_type.eq."mat") then
			read (unitno,NML=matinput)
		endif
		if(TRIM(input_type).eq.'risers')then
			read(unitno,NML=riser)
			!if(volume_fraction.gt.0.1)then
			!   if(I_AM_NODE_ZERO)Write(*,*)'VOLUME FRACTION GREATER THAN 10% FOR RISER FLOWS. STOPPING THE SIMULATION'
			!   PARALLEL_FINISH()
			!   STOP
			!end if
			lybyd = widthbyd
		endif
#if 0
		IF(ISCALON.EQ.1) THEN 
			read(unitno,NML=scal_propt)
		end IF
#endif

		if (iturbon) read(unitno,NML=turb)

		read(unitno,NML=postlist)
		!---------------------------
		close(unitno,STATUS='keep')

		if(my.eq.undefined_I)then
			if(lybyd.eq.undefined_R.or.dbydx.eq.undefined_R)then
				IF (I_AM_NODE_ZERO) PRINT*,'MY UNDEFINED: AND SO IS EITHER LYBYD (',lybyd,') OR DBYDX 	(',DBYDX,'): CHECK THE PARAMETERS'
				IF (I_AM_NODE_ZERO) PRINT*, "STOPPING THE SIMULATION: BBYE U DUD"
				PARALLEL_FINISH()
				STOP
			endif
			!IF (I_AM_NODE_ZERO) PRINT*,'MY UNDEFINED: WILL BE SET IN GENERATE_CONFIGURATION:'
		else
#if 0
			if(mz.eq.undefined_I) then 
				mz = my
				if(I_AM_NODE_ZERO)PRINT*,'MZ UNDEFINED: SETTING EQUAL TO MY, WHICH IS:', MY
			end if

			IF(mx.eq.undefined_I) then 
				mx = my + 1
				IF(I_AM_NODE_ZERO)PRINT*,'MX UNDEFINED: SETTING EQUAL TO MY+1, WHICH IS:', MX
			end IF

			if(I_AM_NODE_ZERO)then
				WRITE(*,*) 'BOX SIZE   = ', mx, my, mz
			end if
#endif
		end if
	END SUBROUTINE read_nmls


	SUBROUTINE write_input(fileno)
		implicit none 
		Integer, Intent(in):: fileno

		call separator(fileno,43,'*')
		write(fileno,*)'INPUT PARAMETERS INCLUDING DEFAULT VALUES'
		call separator(fileno,43,'*')

		call blankline(fileno)
		write(fileno,NML=floandgrid)

		call blankline(fileno)
		write(fileno,NML=part_propt)

		call blankline(fileno)
		write(fileno,NML=Flags)

		call blankline(fileno)
		write(fileno,NML=Files)

		call blankline(fileno)
		write(fileno,NML=gcg)

		call blankline(fileno)
		write(fileno,NML=matinput)

#if 0
		IF(ISCALON.EQ.1) THEN 
			call blankline(fileno)
			write(fileno,NML=scal_propt)
		end IF
#endif
    
	END SUBROUTINE write_input

	SUBROUTINE write_nmls !Just a debug to see how the namelists are
		! written out
		implicit none 
		Logical:: nmlexist
		Integer:: unitno
		CHARACTER*80 :: NODE, FILENAME


		unitno = getnewunit(minunitno, maxunitno)    ! have to use arbitrary unit number

		WRITE(NODE, '(I1)')myid
		FILENAME = TRIM(RUN_NAME)//'_particle_'//TRIM(NODE)//'.opt'

		open(unit=unitno,file=FILENAME,form="formatted",DELIM='APOSTROPHE')
		!Print*,'Writing the namelists to *.opt file'
		!nmlexist = positionnml(unitno,"Properties")

		write(unitno,NML=Floandgrid)

		write(unitno,NML=Part_propt)
		write(unitno,NML=Flags)

		write(unitno,NML=Files)
		!write(unitno,NML=Testing)
		write(unitno,NML=gcg)

		write(unitno,NML=matinput)

		!write(unitno,NML=scal_propt)

		close(unitno, status='keep')
	END SUBROUTINE write_nmls
end module initialize
