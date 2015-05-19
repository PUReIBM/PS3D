PROGRAM mypost_main
  !-----------------------------------------------------------------------
  !	3D code. 
  !	Finite difference in the x-direction
  !	Fourier pseudospectral in the y-direction
  !	Fourier pseudospectral in the z-direction
  ! AUTHOR: RAHUL GARG 
  !	mx	number of interior (u,v,w)-gridpoints in the x-direction
  !	mx1	number of interior (p)-gridpoints in the x-direction!
  !	mxf	number of yz- planes where forcing exists
  !	my2 	number of fourier modes in the y-direction
  !	mz 	number of fourier modes in the z-direction
  !	my	number of (y-) gridpoints on the grid 
  !	mz	number of (z-) gridpoints on the grid 

  !-----------------------------------------------------------------------

#include "../FLO/ibm.h"

	USE precision  
	USE restart_funcs
	USE constants  
	USE global_data
	USE initialize_flo 
!	USE nlcalc
!	USE steptheflow
!	USE initialize  
!	USE dem_mod
!	USE collision_mod
!	USE writeoutput 
	USE machine 
	USE general_funcs
!	USE stepthescalar
!	USE outputscalar
	use mypost_process
	use geometry_2d_3d
	use report
	use fes
	use surface_module
!	use mypost_process_diff_biased
!	use mypost_process_diff_center
!	use physalis_mod
#if PARALLEL
	USE nlarrays, ONLY : ur11, ur22,uatminus1, uatnxp2
	USE nlmainarrays, ONLY : ubcp
#endif
	IMPLICIT NONE

	real(prcn) :: cpu0, cpu1, cputime_used, cputime_hrs, cputime_min
	INTEGER :: run_name_len, orig_rank

	PARALLEL_START()
	GET_NPROCS(comm_group,nproc)

	GET_PROCESSOR_RANK(comm_group,myid)

	! MYID = 0 FOR SERIAL VERSION
	if (I_AM_NODE_ZERO)then
		CALL GETARG(1 , RUN_NAME)
		run_name_len = LEN_TRIM(RUN_NAME)
	end if

	BROADCAST_INT(run_name_len,1,NODE_ZERO,comm_group)
	BROADCAST_CHARARR(RUN_NAME,run_name_len,NODE_ZERO,comm_group)

!	BROADCAST_STRING(RUN_NAME, run_name_len, NODE_ZERO,comm_group)


	CALL calculate_constants

	CALL init_params !set the parameters

	! Virtual processor topology for parallel decomposition
	CREATE_CART_TOPOLOGY(comm_group,1,nproc,xperiodic,.true.,decomp_group)  
	orig_rank = myid

	GET_PROCESSOR_RANK(decomp_group,myid)
	GET_SHIFT_PROCS(decomp_group,0,1,fromproc,toproc)


!	if (irestart/=1) then
!		write (*,*) "IRESTART = ", irestart
!		write (*,*) "CHANGING IT TO 1"
!		irestart = 1
!	endif

	if(I_AM_NODE_ZERO) then
		!--------
		! open error file
		!--------
		eunit = getnewunit(minunitno,maxunitno)
		IF (eunit.LT.0) THEN   ! have to handle this locally
			WRITE(*,*)'Cannot find unused unit for error file'
			WRITE(*,*)'Quitting...'
			PARALLEL_FINISH()
			STOP
		ENDIF
		errfile = TRIM(RUN_NAME)//'_'//TRIM(errfile)
		IF(IRESTART.EQ.0) THEN 
			OPEN(unit=eunit,file=errfile,form="formatted"  &
			   ,status="replace")
		ELSE
			OPEN(unit=eunit,file=errfile,form="formatted"  &
			   ,POSITION="append")
		end IF
		  
		!--------
		! open file for standard out
		!--------

		ounit = getnewunit(minunitno,maxunitno)
		!ounit = -1
		IF (ounit.LT.0) CALL printerror("newunit","ounit")
		!Print*,'UNit for outfile =', ounit, eunit
		outputfile = TRIM(RUN_NAME)//'_'//TRIM(outputfile)
		IF(IRESTART.EQ.0) THEN 
			OPEN(unit=ounit,file=outputfile,form="formatted"  &
				,status="replace")
		ELSE
			OPEN(unit=ounit,file=outputfile,form="formatted"  &
				,POSITION="append")
		end IF
		
		CALL GET_RUN_ID 
		CALL screen_separator(80,'*')

		call separator(ounit,162,'-')
		CALL CPU_TIME (CPU0)
	endif
	BROADCAST_DOUBLE(cpu0,1,NODE_ZERO,decomp_group)

	CALL main
  
CONTAINS 
	SUBROUTINE main
		IMPLICIT NONE  

		INTEGER ::  sflag,intstep,n,i,j,k,ifilstep,idim,lcstart, res_temp, ierr
		INTEGER idum,ifirstcheck,ic,ii, m,snapunit, count_resfiles, iphs
#if PARALLEL
		INTEGER :: nprocold,fnamelen
		LOGICAL :: rstok
#endif
		Integer, SAVE :: unit_t_his

		REAL(prcn) :: dtold,s_temp, re_umax, abtparam, time_since_last_res, CPUTIME_USED_LAST_HRS
		CHARACTER*80 :: conv_file, res_file, filename1
		CHARACTER*50 surface_file
		LOGICAL:: filexist, isopen, CALL_FLOW, stop_criterion

		real(prcn) :: int_dist
		real(prcn) :: taup, tauf, froude

		FROM_POST = .true.
		irestart=1
		if (post_no_flow_mem_alloc) goto 110
		if (isa==1) goto 110
		if (icohesive==1) goto 110

		if(I_AM_NODE_ZERO)then
			unit_t_his  = getnewunit(minunitno,maxunitno)

			conv_file = TRIM(RUN_NAME)//"_CONVERGED"
			res_file = TRIM(RUN_NAME)//"_RESTART"
			!First check if Convergence indicator file exists or not.
			!if it does, then delete it. 
			INQUIRE(FILE=conv_file,EXIST=filexist,OPENED=isopen)
			IF (filexist) THEN
				OPEN(unit = 1000, file=conv_file, status="old")
				close(1000, status="delete")
			end IF
		end if

	
		if(I_AM_NODE_ZERO)then
			!check for this case if something screwed up during the
			! original run and as a result restart files were not properly
			!  written out. Can't happen with the logic in our code. But
			! once hpc4 guys deleted some files during run time  and all this mess
			! hapenned. MIS*_RESTART had 0 written in it and all the
			! restart files were also not written out. RG 11/14/08

			res_temp = 0 
			OPEN(unit = 1000, file=res_file, status="unknown", IOSTAT=IERR)
			READ(1000, *, IOSTAT=IERR) res_temp
#if PARALLEL
			read(1000, *, IOSTAT=IERR) nprocold
			if(nprocold.eq.nproc)rstok = .TRUE.
#endif
			WRITE(*,*)'FIELD IN RESTART INDICATOR FILE =', res_temp, irestart  
			CLOSE(1000, status='keep')

			if(res_temp.eq.0) then 
				CALL screen_separator(80,'E')
				WRITE(*,*) 'WARNING'
				WRITE(*,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
					& 0 or nothing written in it. Something screwed up',' So stopping the post processing'
				WRITE(ounit,*) 'WARNING'
				WRITE(ounit,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
					& 0 or nothing written in it. Something screwed up',' So stopping the post processing'

!				PARALLEL_FINISH()
!				STOP
!#if PARALLEL
!				RSTOK = .TRUE.
!#endif
		   end if

!#if PARALLEL
!			BROADCAST_LOGICAL(rstok,1,NODE_ZERO,decomp_group)
!			BROADCAST_INT(irestart, 1, NODE_ZERO, decomp_group)
!			BROADCAST_INT(iscal_restart, 1, NODE_ZERO, decomp_group)
!			if(.not.rstok)then
!				WRITE(*,'(A,i2)')'RESTARTING RUN WITH DIFFERENT NUMBER OF NODES: ', nproc
!				WRITE(*,'(A,i2,A)')'PREVIOUSLY RUN WITH ', nprocold, ' NODES. RESTART AGAIN.'
!
!				PARALLEL_FINISH()
!				STOP
!			end if
!#endif
		END if


		CALL initflo

		if (POST_NO_FLOW_MEM_ALLOC) goto 110
		if (igeometry==1) goto 110

		call calc_velreal(u, umean, ur)		
		call calc_pressure
		call compute_mean_fluid_velocity
		call fes_compute_new_timestep(.false.)
		call compute_hydrodynamic_forces
		call report_force

110	continue

!		if (I_AM_NODE_ZERO) then
!			call correct_time
!			call combine_history
!		endif
!		call relative_acceleration
		call write_drag
		surface_file = "FES"
		call surface_field(surface_file)
!		call cluster_characteristics
!		call simulated_annealing
!		if (icohesive==1) call cohesive_particles
!		if ((.not.post_no_flow_mem_alloc).and.igeometry==1.and.isa==0) call save_part_restart
!		call fluid_particle_acceleration
!		if (post_no_flow_mem_alloc) call compute_gofr_avg
!		call velocity_output
!		call indicator_output
!		call projection_2d
!		call indicator_3d
!		call calc_part_statistics(1)
!		call reynolds_stress_tensor
!		call cohesive_particles
!		call compute_sijsij
!		call compute_interphase_transfer
!		call c_k
!		call compute_AiVi
!		call pi_groups
!		if (post_no_flow_mem_alloc) call combine_history
!		call output
!		call physalis
!		call post_diff_biased
!		call post_diff_center
!		call dissip_continuous
!		call compute_gofr
!		call uiui_correlation
!		call interstitial_dist(xc,radbdy, int_dist)
		call velocity_output
!		call flow_snapshot2
!		call post_part_stat
!		call volumetric_drag2
!-----------------------------------------


		if(I_AM_NODE_ZERO)then
			WRITE(*,*) 'HERE BECAUSE POST PROCESS ENDED'

			CALL GET_RUN_ID 
			CALL screen_separator(80,'*')
			WRITE(*, '(A40,i2,A,i2,A,i4)') 'SIMULATION ENDED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
			WRITE(*, '(A40,i2,A,i2,A,i2)') 'AT TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND

			call separator(ounit,62,'-')
		end if
		PARALLEL_FINISH()
		STOP
	END SUBROUTINE main
END PROGRAM mypost_main

SUBROUTINE part_snapshot
  USE precision
  USE global_data
  USE general_funcs
  USE dependent_functions

  IMPLICIT NONE
  Integer  :: m
  LOGICAL, SAVE :: first_time_here=.TRUE.
  REAL(prcn) :: ucg
  CHARaCTER*80 :: FILENAME
  CHARACTER(LEN=80) :: formfile
  
  formfile='formatted'
      
  !if(irestart.eq.1)first_time = .FALSE.
  
  IF(first_time_here)THEN
     FILENAME = TRIM(RUN_NAME)//'_part_snapshot.dat'
     CALL  RUN_TIME_FILE_OPENER(partunit,FILENAME, formfile)
     first_time_here = .FALSE.
  END IF
  
  WRITE(partunit,*)'ZONE'
  WRITE(partunit,*)t
  DO m=1,nbody
     WRITE(partunit,'(10(2x,f12.8))')  xc(m,1), xc(m,2), xc(m,3), radbdy(m),velbdy(m,1:ndim),force(m,1:ndim)
  enddo
  
END SUBROUTINE part_snapshot

