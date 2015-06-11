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


PROGRAM dtibm_main
  !-----------------------------------------------------------------------
  !	3D code. 
  !	Fourier pseudospectral in the x-direction
  !	Fourier pseudospectral in the y-direction
  !	Fourier pseudospectral in the z-direction
  !	mx	number of gridpoints in the x-direction
  !	my	number of gridpoints in the y-direction
  !	mz	number of gridpoints in the z-direction
  !-----------------------------------------------------------------------

#include "ibm.h"
	use precision  
	use restart_funcs
	use constants  
	use global_data
	use initialize_flo
	use nlcalc
	use steptheflow
	use initialize  
	use dem_mod
	use collision_mod
	use writeoutput 
	use machine 
	use general_funcs
	!use stepthescalar
	!use outputscalar
	!use fes
	use mypost_process
	use parallel
	implicit none

	real(prcn) :: cpu0, cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cputime_USED, cputime_hrs, cputime_min
	integer :: run_name_len, orig_rank

	!if (I_AM_NODE_ZERO) then
		call GETARG(1, run_name)
		!run_name = "MIS1"
		run_name_len = LEN_trim(run_name)
	!endif

	call calculate_constants
	call init_params !set the parameters

	call initialize_parallel

	! MYID = 0 FOR SERIAL VERSION
	!if (I_AM_NODE_ZERO) then
	!	call GETARG(1, run_name)
	!	run_name_len = LEN_trim(run_name)
	!endif
	BROADCAST_INT(run_name_len, 1, node_zero, comm_cart_2d)
	BROADCAST_CHARARR(run_name, run_name_len, node_zero, comm_cart_2d)

	if (I_AM_NODE_ZERO) then
		!--------
		! open error file
		!--------
		eunit = getnewunit(minunitno,maxunitno)
		if (eunit<0) then   ! have to handle this locally
			write (*,*) 'Cannot find unUSED unit for error file'
			write (*,*) 'Quitting...'
			PARALLEL_FINISH()
			stop
		endif

		errfile = trim(run_name)//'_'//trim(errfile)
		if (IRESTART==0) then 
			open (unit=eunit,file=errfile,form="formatted",status="replace")
		else
			open (unit=eunit,file=errfile,form="formatted",position="append")
		endif
       
		!--------
		! open file for standard out
		!--------
		ounit = getnewunit(minunitno,maxunitno)
		!ounit = -1
		if (ounit<0) call printerror("newunit","ounit")
		!Print*,'UNit for outfile =', ounit, eunit
		outputfile = trim(run_name)//'_'//trim(outputfile)
		if (IRESTART==0) then
			open (unit=ounit,file=outputfile,form="formatted",status="replace")
		else
			open (unit=ounit,file=outputfile,form="formatted",position="append")
		endif
     
		call get_run_id
		call screen_separator(80,'*')
		!call separator(ounit,80,'-')
		if (IRESTART==0) then
#if PARALLEL
			write (*,'(A,i,2A)') 'STARTING A NEW PARALLEL RUN ON ',nproc,' NODES : run_name IS ', run_name
			write (*,'(2(A,i))') 'NPROCY = ', nprocy, 'NPROCZ = ', nprocz
			!write (ounit,'(A,A)') 'STARTING A NEW PARALLEL RUN ON ',nproc,' NODES : run_name IS ', run_name
#else
			write (*,'(A,A)') 'STARTING A NEW RUN: run_name IS ', run_name
			!write (ounit,'(A,A)') 'STARTING A NEW RUN: run_name IS', run_name
#endif
		else
#if PARALLEL
			write (*,'((A,a,a,i,a))') 'RESTARTING AN OLD RUN HAVING RUN_NAME ', run_name,' ON ',nproc,'  NODES'
			write (*,'(a,i,a,i)') 'NPROCY = ', nprocy, 'NPROCZ = ', nprocz
			!write (ounit,'(A,i2,A)') 'RESTARTING AN OLD RUN HAVING run_name = ', run_name,' ON ', nproc, ' NODES'
#else
			write (*,'(A,A)') 'RESTARTING AN OLD RUN HAVING RUN_NAME = ', run_name
			!write (ounit,'(A,A)') 'RESTARTING A OLD RUN HAVING run_name = ', run_name
#endif
		endif
     
		write (*, '(A40,i2,A,i2,A,i4)') 'SIMULATION STARTED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
		write (*, '(A40,i2,A,i2,A,i2)') 'SIMULATION STARTING TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND
		write (*,'(A40,A,A,I1)') 'RUNNING ON MACHINE ', ID_NODE,' WITH RANK ', MYID
		write (*,"(1a)") "USING THE FFTW3-PUReIBM CODE"

		!call separator(ounit,80,'-')

		!write (ounit, '(A40,i2,A,i2,A,i4)') 'SIMULATION STARTED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
		!write (ounit, '(A40,i2,A,i2,A,i2)') 'SIMULATION STARTING TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND
		!write (ounit,'(A40,A)') 'RUNNING ON MACHINE ', ID_NODE

		call screen_separator(80,'*')
		!call separator(ounit,80,'-')
		call CPU_TIME (CPU0) 
	endif
  
	BROADCAST_DOUBLE(cpu0,1,node_zero,comm_cart_2d)
#if PARALLEL
	!write (*,'(3(2x,A25,i4))')'RANK = ', myid, 'SORUCE = ', fromproc, 'DEST = ', toproc
#endif
	BARRIER(comm_cart_2d)
	call main
	!write (*, *) 'Simulation started on date:', date,' at time:', time
  
contains 

	subroutine main
		use init_turb, only : forced_turbulence, statistics
		implicit none  

		integer ::  sflag,intstep,n,i,j,k,res_temp, ierr
		integer idum,ifirstcheck,ic,ii, m,snapunit, count_resfiles, iphs
		integer :: nprocold, nprocyold, nproczold, fnamelen
		logical :: rstok
		integer, save :: unit_t_his

		real(prcn) :: dtold,s_temp, re_umax, abtparam, time_since_last_res, CPUTIME_USED_LAST_HRS
		character*100 :: conv_file, res_file , filename1
		complex(prcn) :: usumloc, usum
		logical:: filexist, isopen, call_flow, killjob, stop_criterion
		real(prcn) ::  rhoofs(ndim),mean_vel(ndim), lag_str_func(ndim), delta_meshpos(ndim)
		real(prcn) :: iterstart_time, iterend_time, avg_iter_time, stat_time
#if PARALLEL
		real(prcn) :: global_avg_iter_time, global_max_iter_time
#endif
		!real(prcn) :: taup, tauf, froude
 
		if (I_AM_NODE_ZERO) then
			unit_t_his  = getnewunit(minunitno,maxunitno)

			conv_file = trim(run_name)//"_CONVERGED"
			res_file = trim(run_name)//"_RESTART"
			!First check if Convergence indicator file exists or not.
			!if it does, then delete it. 
			inquire(file=conv_file, exist=filexist, opened=isopen)
			if (filexist) then
				open (unit = 1000, file=conv_file, status="old")
				close(1000, status="delete")
			endif
		endif


		IF(irestart.eq.0) then 
			if(I_AM_NODE_ZERO)then
				OPEN(unit = 1000, file=res_file, status="unknown")
				res_temp = 0
				write(1000, *) res_temp, nproc, nprocy, nprocz
				close(1000,status="keep")
			end if
		else if(irestart.eq.1)then
			if(I_AM_NODE_ZERO)then
				! check for this case if something screwed up during the
				! original run and as a result restart files were not properly
				! written out. Can't happen with the logic in our code. But
				! once hpc4 guys deleted some files during run time  and all this mess
				! hapenned. MIS*_RESTART had 0 written in it and all the
				! restart files were also not written out. RG 11/14/08

				res_temp = 0 
				OPEN(unit = 1000, file=res_file, status="unknown", IOSTAT=IERR)
				READ(1000, *, IOSTAT=IERR) res_temp, nprocold, nprocyold, nproczold
				if (nprocold==nproc .and. nprocyold==nprocy .and. nproczold==nprocz) rstok = .TRUE.

				WRITE(*,*)'FIELD IN RESTART INDICATOR FILE =', res_temp, irestart  
				CLOSE(1000, status='keep')
				!WRITE(*,*)'IERR =  ', IERR
				if(res_temp.eq.0) then 
					CALL screen_separator(80,'E')
					WRITE(*,*) 'WARNING'
					WRITE(*,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
							& 0 or nothing written in it. Something screwed up',' So starting &
							&the run from scratch'
					WRITE(ounit,*) 'WARNING'
					WRITE(ounit,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
							& 0 written in it. Something screwed up',' So starting &
							&the run from scratch'
							IRESTART = 0
							ISCAL_RESTART = 0
							rstok = .true.

					!Now delete this file and write 0 in it. This is to make sure
					! that in case nothing was writtin in this file, then at the
					! restart stage we do not get any read error. RG 11/16/08

					OPEN(unit = 1000, file=res_file, status="unknown", IOSTAT=IERR)
					CLOSE(1000, status='delete')
					OPEN(unit = 1000, file=res_file, status="unknown")
					write(1000, *) res_temp, nproc, nprocy, nprocz
					close(1000,status="keep")

					CALL screen_separator(80,'E')
				endif
			endIF

			BROADCAST_LOGICAL(rstok,1,NODE_ZERO,decomp_group)
			BROADCAST_INT(irestart, 1, NODE_ZERO, decomp_group)
			BROADCAST_INT(iscal_restart, 1, NODE_ZERO, decomp_group)

			BROADCAST_INT(res_temp, 1, NODE_ZERO, decomp_group)
			BROADCAST_INT(nprocold, 1, NODE_ZERO, decomp_group)
			BROADCAST_INT(nprocyold, 1, NODE_ZERO, decomp_group)
			BROADCAST_INT(nproczold, 1, NODE_ZERO, decomp_group)


			if (.not.rstok) then
				if(I_AM_NODE_ZERO)then
					WRITE(*,'(A,3i4)')'RESTARTING RUN WITH DIFFERENT NUMBER OF NODES: ', nproc, nprocy, nprocz
					WRITE(*,'(A,3i4,A)')'PREVIOUSLY RUN WITH ', nprocold, nprocyold, nproczold, ' NODES. RESTART AGAIN.'
				endif
				PARALLEL_FINISH()
				stop
			endif
		endif


		from_post = .false.
		post_no_flow_mem_alloc = .false.

		call initflo
		if (iturbon) stat_time = eddy_time_i

		final_restart = .false.

		!-------------------------------------------------------------------
		!	timestep through simulation
		!------------------------------------------
		t = 0.d0
		abtparam = 0 
		idumstep = 0
		iglobstep = 0

		t = tstart  !in case of restart, tstart will be equal to the last
		! time step. Restart reading is done from initflo.f90

		if (I_AM_NODE_ZERO) PRINT*,' TOLERANCE FOR FERROR_HIST =',  TOL_FERROR
		call_flow = .TRUE.

		!if (Re==ZERO) then 
		!   ferror_hist = 0.1d0*tol_ferror
		!endif

		flow_converged = ferror_hist.lt.tol_ferror
#if 0
		scal_converged = nu_error_hist.lt.tol_ferror
#endif

		flow_converged = .TRUE.
		!!$grant_converged = .false.
		!!$if (imove==1)grant_converged = .TRUE.

		do iphs = 1, nphases
			flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
		enddo
		!!$if (imove==1) then
		!!$do iphs = 1, nphases
		!!$grant_converged = grant_converged.and.(phase_array(iphs)%gran_error_hist.lt.tol_gran_error)
		!!$enddo
		!!$endif

		call_flow = .not.flow_converged 
		if (iscalon==0) scal_converged = .TRUE.
		if (I_AM_NODE_ZERO) then 
			write (*,*)'FLOW_CONVERGED  = ', flow_CONVERGED
			write (*,*)'SCAL_CONVERGED  = ', SCAL_CONVERGED
			!write (*,*)'GRAN_CONVERGED  = ', GRANT_CONVERGED
			PRINT*,'MAX WALL TIME ALLOWED =', wtime_maxhrs
		endif
		!dO WHILE (t.LE.tend)
		cputime_hrs = zero
		time_since_last_res = zero
		count_resfiles = 0

		!BARRIER(comm_cart_2d)
		stop_criterion = flow_converged.and.scal_converged

		if (imove==1) stop_criterion = .false. ! For a moving particle case run until wall time exceeds the limit.
		if (imove==1) call_flow = .TRUE. ! For a moving particle case always call the flow.

		avg_iter_time = zero
		do while(.not.stop_criterion)
			if (I_AM_NODE_ZERO) call CPU_TIME (CPU2)

			if (I_AM_NODE_ZERO) then
				write (*,'(A,2(2x,g12.5))')'FERROR_HIST AND TOL = ', ferror_hist, TOL_FERROR
				do iphs = 1, nphases
					write (*,'(A,I8,2(2x,g12.5))')'FERROR_HIST AND TOL = ', iphs,phase_array(iphs)%ferror_hist, TOL_FERROR
					if (move_particles) write (*,'(A,I8,2(2x,g12.5))') 'GRAN_ERROR_HIST AND TOL = ', &
												& 		iphs,phase_array(iphs)%gran_error_hist, TOL_GRAN_ERROR
				enddo
			endif
#if 0
			if (ISCALON==1) then
				if (I_AM_NODE_ZERO) write (*,'(A25,g12.5)')'NU ERROR HIST = ', nu_error_hist 
			endif
#endif

			if (I_AM_NODE_ZERO) then
				call screen_separator(80,'-')
				call screen_separator(80,'-')
			endif

			idumstep = idumstep + 1
			iglobstep = iglobstep + 1
       
#if !PARALLEL
			call CPU_TIME (iterstart_time)
#else
			iterstart_time = MPI_WTIME()
#endif
			DO intstep=1,itrmax
				if (imove==1.and.fixed_moving) then
					if ((.not.move_particles))  move_particles = flow_converged
				elseif (imove==1.and.(.not.fixed_moving)) then
					!if ((.not.move_particles))  move_particles = .true.
					move_particles = .true.
				endif

				if (I_AM_NODE_ZERO) write (*,'(A,2(2x,i5))') 'RKS AND GLOBAL STEP # = ', intstep, iglobstep

				if (use_fes) then
					!call fes_solver(intstep)
					!flow_converged = ferror_hist.lt.tol_ferror
					!flow_converged = .TRUE.
					!do iphs = 1, nphases
					!	flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
					!enddo
					!call_flow = .not.flow_converged 
					!if (move_particles) call_flow = .TRUE.
				else
					if (I_AM_NODE_ZERO) write (*,'(A)') 'COMPUTING THE NEW TIME STEP'
					call compute_new_timestep(intstep)

					if (iglobstep==1 .or. mod(iglobstep,skip_num)==0) then
						call reynolds_stress_tensor
						call compute_sijsij
					endif

					if ((mean_vel_to_particles).and.(imove.ne.1).and.(.not.movingcv)) then
						delta_meshpos(1:ndim) = mesh_vel(1:ndim)*dt/dx
						write (*,*)'Delta_meshpos : ', delta_meshpos(1:ndim)
						call interpolate_fields_to_new_mesh(delta_meshpos(1:ndim))
					endif
!!$
!!$          	if (MOD(iglobstep,saveitns)==0) then 
!!$             
!!$            	 call output
!!$             	stop
!!$          		endif

					if (call_flow) then
						call nonlinear(intstep)
						!In non linear now the dt is reset to newly calculated dt based on 
						!cfl criteria
						!Now t is incremented in nl to ensure correct value goes in bcset

						if (I_AM_NODE_ZERO) call CPU_TIME (CPU3)
						if (.not.only_dem) call velstep(sflag,intstep)
						if (I_AM_NODE_ZERO) then
							call CPU_TIME (CPU4)
							vp_time_current = (cpu4-cpu3)
							vp_time = vp_time + vp_time_current

							filename1 = trim(run_name)//'_dt_history'
							unit_t_his=1
							call instant_file_opener(filename1,unit_t_his,.true.)
							write (unit_t_his, '(20(2x,g17.8))') real(iglobstep,prcn), DT, (tend-t)/dt + real(iglobstep, PRCN), uchar(1)*dt/dx
							close(unit_t_his)
						endif

						!if (debug_check) then
						!	usumloc = SUM(u(1:nx,1,1,1))
						!	GLOBAL_complex_SUM(usumloc,usum,1,comm_cart_2d)
						!	if (I_AM_NODE_ZERO) then
						!		write (*,'(A25,2(2x,g12.5))') 'UAVG = ', usum/mx1 !SUM(u(1:mx1,1,1,1))/mx1 
						!		write (*,'(A25,2(2x,g12.5))') 'MAX DIVERGENCE = ', divmax
						!	endif
						!endif
						!if (t.gt.tendUSED)
						flow_converged = ferror_hist.lt.tol_ferror
						flow_converged = .TRUE.
						do iphs = 1, nphases
							flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
						enddo
					    
						call_flow = .not.flow_converged 
						if (move_particles) call_flow = .TRUE.
					else
						if (I_AM_NODE_ZERO) write (*,*)'NOT CALLING FLOW ROUTINES BECAUSE EITHER call_flow = ', &
																&	call_flow, 'OR flow_CONVERGED = ', flow_converged
					endif

#if 0
					if (iscalon==1) then
						if (I_AM_NODE_ZERO) call CPU_TIME (CPU5)
#if PARALLEL
						if (flow_converged) then
							if (.not.vel_converted) then
								call ff2cr(u(2,1:my2,1:mz,1),ur11(1:my,1:mz))
								call ff2cr(u(nx-1,1:my2,1:mz,1),ur22(1:my,1:mz))
								RSENDRECV(ur11(1,1),my*mz,fromproc,1,uatnxp2(1,1),my*mz,toproc,1,comm_cart_2d,status)
								RSENDRECV(ur22(1,1),my*mz,toproc,1,uatminus1(1,1),my*mz,fromproc,1,comm_cart_2d,status)
								uatnxp2(:,:) = uatnxp2(:,:)+umean(1)
								uatminus1(:,:) = uatminus1(:,:)+umean(1)
								vel_converted = .TRUE.
							endif
						else
							uatminus1(1:my,1:mz) = ubcp(-1,1:my,1:mz,1)
							uatnxp2(1:my,1:mz) = ubcp(nx+2,1:my,1:mz,1)
						endif
#endif
						!nlphi now called in scalstep 
						call scalstep(sflag,intstep)

						!if (t.gt.tendUSED)
						scal_converged = nu_error_hist.lt.tol_ferror
						!scal_converged  = .true.
						if (I_AM_NODE_ZERO) call CPU_TIME (CPU6)
					endif
#endif
				endif

				stop_criterion = flow_converged !.and.scal_converged
				if (imove==1) stop_criterion = .false.
#if 0
				if (ISCALON==1) then
					usumloc = SUM(phif (1:nx,1,1,1))
					GLOBAL_complex_SUM(usumloc,usum,1,comm_cart_2d)
					if (I_AM_NODE_ZERO) write (*,'(A25,2(2x,g12.5))') 'PHI_AVG = ', usum/mx1
				endif
#endif
				if (I_AM_NODE_ZERO)call screen_separator(80,'*')
				!first_pass = .false.   
				!cf = cforig
			enddo

			if (mod(iglobstep,saveitns)==1) then
				!call output
				!if (iscalon==1) call output_scal
				if (write_snapshot) call flow_snapshot
				!call compute_AiVi
			endif

			if (I_AM_NODE_ZERO) then
				call CPU_TIME (CPU1)
				!BROADCAST_DOUBLE(cpu1,1,node_zero,comm_cart_2d)

				CPUTIME_USED_LAST_HRS = CPUTIME_HRS
				CPUTIME_USED  = CPU1 - CPU0
				CPUTIME_HRS = CPUTIME_USED/(3600.d0)
				CPUTIME_MIN = CPUTIME_USED/(60.d0)
				time_since_last_res = time_since_last_res + (-cputime_USED_last_hrs+cputime_hrs)

				PRINT*,'CPU TIMEHRS = ', CPUTIME_HRS

				hydro_time_current = (cpu1-cpu2)-dem_time_current
				hydro_time = hydro_time + hydro_time_current

				scalar_time_current = (cpu6-cpu5)
				scalar_time = scalar_time + scalar_time_current

				write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN HYDRO   (S), CURRENT, TOTAL, AVG : ', hydro_time_current , hydro_time, hydro_time/iglobstep
				if (use_fes) write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN FOURIER(S), CURRENT, TOTAL, AVG : ', fourier_time_current , fourier_time, fourier_time/iglobstep
				write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN TIMESTEP(S), CURRENT, TOTAL, AVG : ', new_timestep_time_current , new_timestep_time, new_timestep_time/iglobstep
				write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN NL      (S), CURRENT, TOTAL, AVG : ', nl_time_current , nl_time, nl_time/iglobstep
				write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN BC      (S), CURRENT, TOTAL, AVG : ', bc_time_current , bc_time, bc_time/iglobstep
				!write (*,*) "********************"
				!write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN MPG      (S),                    : ', mpg_time
				!write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN PGRAD    (S),                    : ', pgrad_time
				!write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN REVERSAL (S),                    : ', reversal_time
				!write (*,*) "********************"
				write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN VP      (S), CURRENT, TOTAL, AVG : ', vp_time_current , vp_time, vp_time/iglobstep
				if (use_fes) then
					hydroforce_time = hydroforce_time + hydroforce_time_current
					write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN HF (S)    , CURRENT, TOTAL, AVG : ', &
							& hydroforce_time_current , hydroforce_time, hydroforce_time/iglobstep
				endif

				if (move_particles) then
					write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN DEM (S),   CURRENT, TOTAL, AVG : ', dem_time_current, dem_time, dem_time/iglobstep
					write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN H+DEM(S),  CURRENT, TOTAL, AVG : ', dem_time_current+hydro_time_current, &
																										& dem_time+hydro_time, (dem_time+hydro_time)/iglobstep
				endif

				if (iscalon==1) then
					write (*,*)
					write (*,'(A,3(2x,g15.8))') 'CPUTIME USED IN SC (S),    CURRENT, TOTAL, AVG : ', scalar_time_current, scalar_time, &
																																& scalar_time/iglobstep
					write (*,'(A,3(2x,g15.8))') 'CPUTIME USED H+D+S (S),    CURRENT, TOTAL, AVG : ', dem_time_current+hydro_time_current+scalar_time_current, dem_time+hydro_time, (dem_time+hydro_time+scalar_time)/iglobstep
				endif

				write (*,*)
				write (*,'(A,4(2x,g15.8))') 'CPUTIME USED (H:M:S)', CPUTIME_HRS,CPUTIME_MIN, CPUTIME_USED
				write (*,'(A,4(2x,g15.8))') 'TIME SINCE LAST RES, RES_TIME', time_since_last_res, saveforres_time

				if (CPUTIME_HRS.GT.WTIME_MAXHRS) then 
					PRINT*,'KIILING THE JOB BECAuse WTIME_MAXHRS EXCEEDED'
					PRINT*, 'CONVERGENCE NOT YET REACHED'
					killjob = .TRUE.
				else
					killjob = .false.
				endif
			endif

			if (trim(input_type).eq."single-phase") then
				if (iturbon.and.forced_turbulence) then
					if (t>stat_time) then
						stat_time = t+eddy_time_i
						call statistics(1)
					endif
				endif
				if (t>100*eddy_time_i) killjob = .true.
			endif


			BROADCAST_LOGICAL(killjob,1,node_zero,comm_cart_2d)
			BROADCAST_DOUBLE(time_since_last_res,1,node_zero,comm_cart_2d)

!if (iglobstep==5) killjob=.true.
			if (killjob) goto 999

#if !PARALLEL
			call CPU_TIME (iterend_time)
			avg_iter_time = avg_iter_time + (iterend_time - iterstart_time)
#else
			iterend_time = MPI_WTIME()
			avg_iter_time = (iterend_time - iterstart_time)
			GLOBAL_DOUBLE_MAX(avg_iter_time,global_max_iter_time,1,comm_cart_2d)
			global_avg_iter_time = global_avg_iter_time + global_max_iter_time
#endif
			!taup = rhos/rhof * dia_phys**2/18/vis
			!tauf = dia_phys / sqrt(dot_product(ufmean(:)-usmean(:), ufmean(:)-usmean(:)))
			!froude = mpg(1) * taup**2 / doml(1)
			!if (I_AM_NODE_ZERO) then
			!	write (*,"(1a,3d15.7)") "UFMEAN = ", ufmean
			!	write (*,"(1a,3d15.7)") "USMEAN = ", usmean
			!	write (*,"(1a,3d15.7)") "MPG    = ", mpg
			!	write (*,"(1a,1d15.7)") "RE     = ", RE
			!	write (*,"(1a,3d15.7)") "DOML   = ", doml
			!	write (*,"(1a,1d15.7)") "FROUDE = ", froude
			!	write (*,"(1a,3d15.7)") "TAUP, TAUF, ST = ", taup, tauf, taup/tauf
			!endif
			!PARALLEL_FINISH()
			!stop

			!output intermediate restart file for every 100 steps
!!$

			if ((saveforrestart==1).AND.(time_since_last_res.gt.saveforres_time)) then 
				count_resfiles = count_resfiles+1
				if (I_AM_NODE_ZERO) then
					write (*,'(A,g17.8,A,g17.8,A)')'WRITING RESTART fileS BECAuse time_since_last_res (= ', &
						&	 time_since_last_res,') is > saveforres_time (= ', saveforres_time, ')'
					write (ounit,'(A,g17.8,A,g17.8,A)')'WRITING RESTART fileS BECAuse time_since_last_res (= ', &
						&	time_since_last_res,') is > saveforres_time (= ', saveforres_time, ')'
					write (*,'(A,i4)')'NUMBER OF TIMES RESTART fileS WRITTEN IN THIS RUN = ', count_resfiles
					write (ounit,'(A,i4)')'NUMBER OF TIMES RESTART fileS WRITTEN IN THIS RUN = ', count_resfiles
				endif
				call save_part_restart           
				time_since_last_res = zero
			endif
			first_pass = .false.   

		enddo

		if (I_AM_NODE_ZERO) then 
			write (*,*) 'HERE BECAuse CONVERGENCE REACHED'
			write (*,'(A,2(2x,g12.5))')'FERROR_HIST AND TOL = ', ferror_hist, TOL_FERROR
#if 0
			if (ISCALON==1) write (*,'(A25,g12.5)')'NU ERROR HIST = ', nu_error_hist 
#endif
			open (unit = 1000, file=conv_file, status="replace")
			close(1000, status="keep")
		endif

999	continue 

		if (I_AM_NODE_ZERO) then
			open (unit=1, file=trim(run_name)//"_time.dat", status="replace", action="write")
			write (1,"(4i,10d15.7)") nproc, nprocy, nprocz, my, &
			&			hydro_time_current, new_timestep_time_current, nl_time_current, bc_time_current, vp_time_current , &
			&			hydro_time/iglobstep, new_timestep_time/iglobstep, nl_time/iglobstep, bc_time/iglobstep, vp_time/iglobstep
			close (1)
		endif

		if (saveforrestart==1) then 
			count_resfiles = count_resfiles+1
			if (I_AM_NODE_ZERO) then
				write (*,'(A,g17.8)')'WRITING RESTART fileS AT THE END: TIME_SINCE_LAST_RESTART file WRITTEN = ', time_since_last_res
				write (ounit,'(A,g17.8)')'WRITING RESTART fileS AT THE END: TIME_SINCE_LAST_RESTART file WRITTEN = ', time_since_last_res

				write (*,'(A,i4)')'NUMBER OF TIMES RESTART fileS WRITTEN IN THIS RUN = ', count_resfiles
				write (ounit,'(A,i4)')'NUMBER OF TIMES RESTART fileS WRITTEN IN THIS RUN = ', count_resfiles
			endif
			final_restart = .false.
			call save_part_restart
		endif

		!call output



#if 0 
		if (iscalon==1) call output_scal
#endif

111	FORMAT(i4)
		if (I_AM_NODE_ZERO) then
			PRINT*,'Done with writing the output'
			close(unit_t_his,status='keep')
			close(unitnormdrag,status='keep')
			close(unitnormdragchem,status='keep')
			close(unitdragtavg,status='keep')
			close(unitdrag_comps,status='keep')
			close(sphrunit,status='keep')
			if (move_particles) then
				close(unitpartinfo,status='keep')
			endif
			!if (base_state_saved) close(partunit, status='keep')
			!t = 20.0
			!Print*,'time before entering test flux = ', t
			!call calc_anal

			call GET_RUN_ID 
			call screen_separator(80,'*')
			write (*, '(A40,i2,A,i2,A,i4)') 'SIMULATION ENDED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
			write (*, '(A40,i2,A,i2,A,i2)') 'AT TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND

			call separator(ounit,62,'-')

			write (ounit, '(A40,i2,A,i2,A,i4)') 'SIMULATION ENDED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
			write (ounit, '(A40,i2,A,i2,A,i2)') 'AT TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND
#if !PARALLEL       
			write (*, '(A40,2x,g17.5, 2x,A10 )') 'COST/GRIDCELL/ITER:', avg_iter_time/real(iglobstep, prcn), 'SECONDS'
#else
			write (*, '(A40,2x,g17.5, 2x,A10 )') 'COST/GRIDCELL/ITER:', global_avg_iter_time/real(iglobstep, prcn), 'SECONDS'
#endif
			call screen_separator(80,'*')
			call separator(ounit,62,'-')
		endif


		PARALLEL_FINISH()
		stop
	END subroutine main
END PROGRAM dtibm_main

subroutine part_snapshot
	use global_data
	use general_funcs
	use dependent_functions

	implicit none
	integer  :: m
	logical, save :: first_time_here=.TRUE.
	real(prcn) :: ucg
	character*100 :: filename, formfile
	formfile='formatted'

	!if (irestart==1)first_time = .false.
	if (first_time_here) then
		filename = trim(run_name)//'_part_snapshot.dat'
		call  RUN_TIME_file_OPENER(partunit,filename, formfile)
		first_time_here = .false.
	endif

	write (partunit,*)'ZONE'
	write (partunit,*)t
	DO m=1,nbody
		write (partunit,'(10(2x,f12.8))')  xc(m,1), xc(m,2), xc(m,3), radbdy(m),velbdy(m,1:ndim),force(m,1:ndim)
	enddo
END subroutine part_snapshot

#if 0
subroutine u_periodic_bc
	use precision  
	use constants  
	use global_data
	use scalar_data
	implicit none  
	integer :: j, k 
	real(prcn) :: a(4), b(4), c(4), r1(4), sol(4)
	DO k = 1, mz
		DO j = 1, my2 
			uin(j,k,:)  = u(mx1, j,k,:)
			uout(j,k,:)  = u(2, j,k,:)
			if (iscalon==1) then 
				phiin(j,k,:) = phif (mx1,j,k,:)
				phiout(j,k,:) = phif (2,j,k,:)
			endif
		enddo
	enddo
	a(1:4) = one
	c(1:4) = one
	if (I_AM_NODE_ZERO) then
		b(1) = 3.d0
		b(2) = 3.d0
		r1(1) = 1.d0
		r1(2) = 1.d0
	else
		b(1) = 2.d0
		b(2) = 6.d0 
		r1(1) = 3.d0
		r1(2) = 0.5d0 
	endif
!    call tridag(a,b,c,r,sol,4)
!    call mpi_tridag(a(1:2),b(1:2),c(1:2),r1(1:2),sol,2)
END subroutine u_periodic_bc
#endif
