PROGRAM dtibm_main
  !-----------------------------------------------------------------------
  !	3D code. 
  !	Fourier pseudospectral in the x-direction
  !	Fourier pseudospectral in the y-direction
  !	Fourier pseudospectral in the z-direction
  ! AUTHOR: RAHUL GARG, SUDHEET TENNETI, MOHAMMAD MEHRABADI
  !	mx	number of gridpoints in the x-direction
  !	my	number of gridpoints in the y-direction
  !	mz	number of gridpoints in the z-direction
  !-----------------------------------------------------------------------

#include "../FLO/ibm.h"
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
 

		IF(irestart.eq.0) then 
			if(I_AM_NODE_ZERO)then
				write (*,*) "SETTING IRESTART=1..."
			end if
		endif

		if(irestart.eq.1)then
			if(I_AM_NODE_ZERO)then
				!check for this case if something screwed up during the
				! original run and as a result restart files were not properly
				!  written out. Can't happen with the logic in our code. But
				! once hpc4 guys deleted some files during run time  and all this mess
				! hapenned. MIS*_RESTART had 0 written in it and all the
				! restart files were also not written out. RG 11/14/08

				res_file = trim(run_name)//"_RESTART"

				res_temp = 0 
				OPEN(unit = 1000, file=res_file, status="unknown", IOSTAT=IERR)
				READ(1000, *, IOSTAT=IERR) res_temp, nprocold, nprocyold, nproczold
#if PARALLEL
				if (nprocold==nproc .and. nprocyold==nprocy .and. nproczold==nprocz) rstok = .TRUE.
#else
				rstok = .TRUE.
#endif

				WRITE(*,*)'FIELD IN RESTART INDICATOR FILE =', res_temp, irestart  
				CLOSE(1000, status='keep')
				!WRITE(*,*)'IERR =  ', IERR
				if(res_temp.eq.0) then 
					CALL screen_separator(80,'E')
					WRITE(*,*) 'WARNING'
					WRITE(*,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
							& 0 or nothing written in it. Something screwed up',' So starting &
							&the run from scratch'
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


			if (res_temp==0) then
				if(I_AM_NODE_ZERO)then
					WRITE (*,'(A)') "RESTART COUNT = 0. SOMETHING IS WRONG IN RESTART FILE."
					WRITE (*,'(A)') "TERMINATING THE SIMULATION."
				endif
				PARALLEL_FINISH()
				stop
			endif


			if (.not.rstok) then
				if(I_AM_NODE_ZERO)then
					WRITE(*,'(A,3i4)')'RESTARTING RUN WITH DIFFERENT NUMBER OF NODES: ', nproc, nprocy, nprocz
					WRITE(*,'(A,3i4,A)')'PREVIOUSLY RUN WITH ', nprocold, nprocyold, nproczold, ' NODES. RESTART AGAIN.'
				endif
				PARALLEL_FINISH()
				stop
			endif
		endif


		FROM_POST = .true.
		if (post_no_flow_mem_alloc) goto 110
		!if (isa==1) goto 110
		!if (icohesive==1) goto 110

		call initflo
		if (iturbon) stat_time = eddy_time_i

		final_restart = .false.

		!if (POST_NO_FLOW_MEM_ALLOC) goto 110
		if (igeometry==1) goto 110



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

		!flow_converged = ferror_hist.lt.tol_ferror
		!flow_converged = .TRUE.
		!do iphs = 1, nphases
		!	flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
		!enddo
		!call_flow = .not.flow_converged 
		!if (iscalon==0) scal_converged = .TRUE.
		if (I_AM_NODE_ZERO) then 
			write (*,*)'FLOW_CONVERGED  = ', flow_CONVERGED
			write (*,*)'SCAL_CONVERGED  = ', SCAL_CONVERGED
			!write (*,*)'GRAN_CONVERGED  = ', GRANT_CONVERGED
			PRINT*,'MAX WALL TIME ALLOWED =', wtime_maxhrs
		endif

		cputime_hrs = zero
		time_since_last_res = zero
		count_resfiles = 0

		BARRIER(comm_cart_2d)
		stop_criterion = flow_converged.and.scal_converged

		if (imove==1) stop_criterion = .false. ! For a moving particle case run until wall time exceeds the limit.
		if (imove==1) call_flow = .TRUE. ! For a moving particle case always call the flow.

		avg_iter_time = zero


		if (I_AM_NODE_ZERO) call CPU_TIME (CPU2)

		if (I_AM_NODE_ZERO) then
			write (*,'(A,2(2x,g12.5))')'FERROR_HIST AND TOL = ', ferror_hist, TOL_FERROR
			do iphs = 1, nphases
				write (*,'(A,I8,2(2x,g12.5))')'FERROR_HIST AND TOL = ', iphs,phase_array(iphs)%ferror_hist, TOL_FERROR
				if (move_particles) write (*,'(A,I8,2(2x,g12.5))') 'GRAN_ERROR_HIST AND TOL = ', &
											& 		iphs,phase_array(iphs)%gran_error_hist, TOL_GRAN_ERROR
			enddo
		endif

		if (I_AM_NODE_ZERO) then
			call screen_separator(80,'-')
			call screen_separator(80,'-')
		endif
    
#if !PARALLEL
		call CPU_TIME (iterstart_time)
#else
		iterstart_time = MPI_WTIME()
#endif

		if (imove==1.and.fixed_moving) then
			if ((.not.move_particles))  move_particles = flow_converged
		elseif (imove==1.and.(.not.fixed_moving)) then
			!if ((.not.move_particles))  move_particles = .true.
			move_particles = .true.
		endif

		intstep = 1
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

			if ((mean_vel_to_particles).and.(imove.ne.1).and.(.not.movingcv)) then
				delta_meshpos(1:ndim) = mesh_vel(1:ndim)*dt/dx
				write (*,*)'Delta_meshpos : ', delta_meshpos(1:ndim)
				call interpolate_fields_to_new_mesh(delta_meshpos(1:ndim))
			endif

			if (call_flow) then
				call nonlinear(intstep)

				if (I_AM_NODE_ZERO) call CPU_TIME (CPU3)
				!if (.not.only_dem) call velstep(sflag,intstep)
				!if (I_AM_NODE_ZERO) then
				!	call CPU_TIME (CPU4)
				!	vp_time_current = (cpu4-cpu3)
				!	vp_time = vp_time + vp_time_current

				!	filename1 = trim(run_name)//'_dt_history'
				!	unit_t_his=1
				!	call instant_file_opener(filename1,unit_t_his,.true.)
				!	write (unit_t_his, '(20(2x,g17.8))') real(iglobstep,prcn), DT, (tend-t)/dt + real(iglobstep, PRCN), uchar(1)*dt/dx
				!	close(unit_t_his)
				!endif

				!if (debug_check) then
				!	usumloc = SUM(u(1:nx,1,1,1))
				!	GLOBAL_complex_SUM(usumloc,usum,1,comm_cart_2d)
				!	if (I_AM_NODE_ZERO) then
				!		write (*,'(A25,2(2x,g12.5))') 'UAVG = ', usum/mx1 !SUM(u(1:mx1,1,1,1))/mx1 
				!		write (*,'(A25,2(2x,g12.5))') 'MAX DIVERGENCE = ', divmax
				!	endif
				!endif
				!if (t.gt.tendUSED)
				!flow_converged = ferror_hist.lt.tol_ferror
				!flow_converged = .TRUE.
				!do iphs = 1, nphases
				!	flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
				!enddo
			    
				!call_flow = .not.flow_converged 
				!if (move_particles) call_flow = .TRUE.
			else
				if (I_AM_NODE_ZERO) write (*,*)'NOT CALLING FLOW ROUTINES BECAUSE EITHER call_flow = ', &
														&	call_flow, 'OR flow_CONVERGED = ', flow_converged
			endif
		endif



110	continue

		!if (post_no_flow_mem_alloc) call combine_history
		!call post_part_stat
		!call velocity_output
		call write_drag
		call flow_snapshot
		!call pressure_velocity_vs_theta
		!call kp_source_dissipation
		!call fluid_acceleration_pdf
		!call compute_energy_spectrum_function
		!call Apr_App_scatter
		!call reynolds_stress_tensor
		!call compute_sijsij



		if (I_AM_NODE_ZERO) call screen_separator(80,'*')

		if (I_AM_NODE_ZERO) then 
			write (*,*) 'HERE BECASUE POST PROCESS ENDED'
		endif

		if (I_AM_NODE_ZERO) then
			call GET_RUN_ID 
			call screen_separator(80,'*')
			write (*, '(A40,i2,A,i2,A,i4)') 'SIMULATION ENDED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
			write (*, '(A40,i2,A,i2,A,i2)') 'AT TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND

			call screen_separator(80,'*')
		endif


		PARALLEL_FINISH()
		stop
	END subroutine main
END PROGRAM dtibm_main

