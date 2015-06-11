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

MODULE restart_funcs
#include "ibm.h"
	USE precision          ! Independent modules
	USE general_funcs
	USE global_data        ! Dependent modules
	!USE scalar_data 
	USE errormesgs
	Implicit none

	Private
	Public:: save_part_restart, read_restart, delete_file, filename_gen
  
  !-------
CONTAINS
   !-------
  
  !-------
  !-------
  ! Save restart file
  !-------
  !-------
	SUBROUTINE save_part_restart
		IMPLICIT NONE

		INTEGER  :: count_routine
		CHARACTER*10 :: filename1
		CHARACTER*100 :: filename,stat
		LOGICAL :: filexist, isopen
		INTEGER :: runit,strlen,node
		!CALL hello
		!PARALLEL_FINISH()
		!STOP

		If (trim(rstsave).eq."formatted") then
			!call save_formatted
		elseif (trim(rstsave).eq."unformatted") then
			runit = getnewunit(minunitno,maxunitno)
			!BROADCAST_INT(runit,1,NODE_ZERO,comm_cart_2d)
			if (runit.lt.0) call printerror("newunit","runit")
			if(I_AM_NODE_ZERO)then
				filename = TRIM(RUN_NAME)//'_RESTART'
				open(unit=runit,file=trim(filename),form="formatted",status="unknown")    
				read(runit, *) count_routine, nproc, nprocy, nprocz
				close(runit, status='keep')
			endif
			BROADCAST_INT(count_routine,1,NODE_ZERO,comm_cart_2d)
			count_routine = count_routine+1
			IF (count_routine.eq.3) count_routine = 1

			if(I_AM_NODE_ZERO)then
				WRITE(*,*) 'COUNT FOR DUMPING RESTART FILES ', COUNT_ROUTINE
				WRITE(ounit,*) 'COUNT FOR DUMPING RESTART FILES ', COUNT_ROUTINE
			end if

			if (from_post.and.igeometry==1) return

			!CALL filename_GEN(filename,'u',count_routine)
			!INQUIRE(FILE=TRIM(filename),EXIST=filexist,OPENED=isopen)

			!IF (.NOT.filexist) THEN
			!	stat="new"
			!	call save_unformatted(count_routine,stat)
			!elseif(filexist.AND..NOT.isopen) THEN
				stat="replace"
				call save_unformatted(count_routine,stat)
			!ENDIF
		Endif

		count_restart = count_routine

		BARRIER(comm_cart_2d)
		!-------
		! write to output unit
		!-------
   END SUBROUTINE save_part_restart
  
	SUBROUTINE filename_GEN(filename,VAR_NAME,RES_COUNT)
		IMPLICIT NONE

		Character(LEN=*),Intent(out) :: filename
		Character(LEN=*),Intent(in) :: var_name
		INTEGER, INTENT(in) :: RES_COUNT

		CHARACTER*10 :: filename1
		Character*100 :: filenameLOC
		Integer :: node, strlen

		BARRIER(comm_cart_2d)

		if(I_AM_NODE_ZERO)then
#if PARALLEL
			!^^^^ MODIFIED FOR MORE PROCS
			if (nproc<=100) then
				write (filename1,fmt="('NODE',i2.2,'_',i1)") myid, res_count
			elseif (nproc<=1000) then
				write (filename1,fmt="('NODE',i3.3,'_',i1)") myid, res_count
			elseif (nproc<=10000) then
				write (filename1,fmt="('NODE',i4.4,'_',i1)") myid, res_count
			elseif (nproc<=100000) then
				write (filename1,fmt="('NODE',i5.5,'_',i1)") myid, res_count
			elseif (nproc<=1000000) then
				write (filename1,fmt="('NODE',i6.6,'_',i1)") myid, res_count
			endif
#else
			WRITE(filename1, '(I1)') res_count
#endif
			filename = TRIM(RUN_NAME)//'_'//TRIM(var_name)//'_'//TRIM(filename1)//'.rst'
			filenameLOC = ""
		else
			filename = ""
		end if

		if(I_AM_NODE_ZERO)then
			do node=1,nproc-1
				!^^^^ MODIFIED FOR MORE PROCS
				if (nproc<=100) then
					write (filename1,fmt="('NODE',i2.2,'_',i1)") node, res_count
				elseif (nproc<=1000) then
					write (filename1,fmt="('NODE',i3.3,'_',i1)") node, res_count
				elseif (nproc<=10000) then
					write (filename1,fmt="('NODE',i4.4,'_',i1)") node, res_count
				elseif (nproc<=100000) then
					write (filename1,fmt="('NODE',i5.5,'_',i1)") node, res_count
				elseif (nproc<=1000000) then
					write (filename1,fmt="('NODE',i6.6,'_',i1)") node, res_count
				endif

				!filenameLOC = TRIM(RUN_NAME)//'_'//TRIM(var_name)//'_'//TRIM(filename1)//'.rst'
				!!PRINT*,'NODE ZERO', node, filenameloc
				!SEND_STRING(filenameloc,strlen,node,0,1,comm_cart_2d)


				filenameLOC = TRIM(RUN_NAME)//'_'//TRIM(var_name)//'_'//TRIM(filename1)//'.rst'
				strlen = LEN_trim(filenameLOC)

				SEND_INT(strlen,1,node,0,comm_cart_2d)
				SEND_CHARACTER(filenameloc,strlen,node,1,comm_cart_2d)
			enddo
		else  
			!RECV_STRING(filename,strlen,node_zero,0,1,comm_cart_2d,status)
			RECV_INT(strlen,1,node_zero,0,comm_cart_2d,status)
			RECV_CHARACTER(filename,strlen,node_zero,1,comm_cart_2d,status)
		endif

		BARRIER(comm_cart_2d)
	END SUBROUTINE filename_GEN


  !-------
  !-------
  ! If restart file is unformatted
  !-------
  !-------

	SUBROUTINE save_unformatted(count_routine,stat)
		use fftw3_interface
		Implicit none

		INTEGER, Intent(in) :: count_routine

		CHARACTER*100,Intent(in) :: stat
		INTEGER :: count_tmp, strlen, node, iphs, nerr_steps_tmp
		CHARACTER*100 :: filename
		CHARACTER*100 :: filename_config, filename_neighb, filename_u, filename_p, filename_nl, filename_scal
		CHARACTER*10 :: filename1
		LOGICAL:: filexist, isopen


		runit = getnewunit(minunitno,maxunitno)
		if (runit.lt.0) call printerror("newunit","runit")
		if(I_AM_NODE_ZERO)then
			WRITE (filename1, '(I1)') count_routine 	
		endif

		!-------
		! Basic grid data
		!-------
		!    PRINT*,'sat = ', myid, stat
		if (I_AM_NODE_ZERO) then
			filename_config = TRIM(RUN_NAME)//'_sphere_config.rst'
			open(unit=runit,file=trim(filename_config),form=trim(rstsave),status="unknown")

			write(runit) nbody
			if (nbody>0) then
				write(runit)nphases
				do iphs = 1, nphases
					write(runit) phase_array(iphs)%npart
				end do
				do iphs = 1, nphases
					write(runit) phase_array(iphs)%dia
				end do
			endif

			write(runit) DOML(1:ndim)
			write(runit) mx,my,mz
			write(runit) dx,dy,dz

			if (nbody>0) write(runit) move_particles

			write(runit) char_length

			if (nbody>0) then
				write(runit) xc(1:nbody, 1:3)
				write(runit) radbdy(1:nbody)
				write(runit) velbdy(1:nbody, 1:3)
				write(runit) frame_vel(1:ndim)
				write(runit) frame_pos(1:ndim)
			endif
			close(runit, status='keep')

			if (imove.eq.1) then
				filename_config = TRIM(RUN_NAME)//'_sphere_config_'//TRIM(filename1)//'.rst'
				open (unit=runit,file=trim(filename_config),form=trim(rstsave),status="unknown")
				write (runit) nbody
				if (nbody>0) then
					write (runit) nphases
					do iphs = 1, nphases
						write(runit) phase_array(iphs)%npart
					end do
					do iphs = 1, nphases
						write(runit) phase_array(iphs)%dia
					end do
				endif

				write(runit) DOML(1:ndim)
				write(runit) mx,my,mz
				write(runit) dx,dy,dz

				if (nbody>0) write(runit) move_particles

				write(runit) char_length

				if (nbody>0) then
					write(runit) xc(1:nbody, 1:ndim)
					write(runit) radbdy(1:nbody)
					write(runit) velbdy(1:nbody, 1:ndim)
					write(runit) frame_vel(1:ndim)
					write(runit) frame_pos(1:ndim)
				endif
				close(runit, status='keep')
			endif
		endif

#if PARALLEL
		CALL filename_GEN(filename_neighb,'neighb',count_routine)
		open(unit=runit,file=trim(filename_neighb),form="formatted",status=stat) 
		write(runit,"(5i)") myid, fromy, toy, fromz, toz
		write(runit,"(3i)") local_ni(1), local_ni(2), local_ni(3)
		write(runit,"(3i)") local_no(1), local_no(2), local_no(3)
		close(runit, status='keep') 
#endif

		CALL filename_GEN(filename_u,'u',count_routine)
		open(unit=runit,file=trim(filename_u),form=trim(rstsave),status=stat)

		write(runit) t, dt
		write(runit) umean(1), umean(2), umean(3)
		write(runit) usmean_des(1), usmean_des(2), usmean_des(3)
		write(runit) ufmean_des(1), ufmean_des(2), ufmean_des(3)
		write(runit) cf 
		write(runit) nerr_steps
		write(runit) ferror, fold
		write(runit) ferror_array(1:nerr_steps)
		do iphs = 1, nphases
			write(runit) phase_array(iphs)%ferror, phase_array(iphs)%fold
			if (imove.eq.1) write(runit) phase_array(iphs)%grant_error, phase_array(iphs)%grant_old
		end do

		do iphs = 1, nphases
			write (runit) phase_array(iphs)%ferror_array(1:nerr_steps)
			if (imove.eq.1) write(runit) phase_array(iphs)%grant_array(1:nerr_steps)
		end do

		write (runit) u
		close (runit, status='keep')

		CALL filename_GEN(filename_nl,'nl',count_routine)
		open (unit=runit,file=trim(filename_nl),form=trim(rstsave),status=stat)
		write (runit) nl
		close (runit, status='keep')

		CALL filename_GEN(filename_p,'p',count_routine)
		open (unit=runit,file=TRIM(filename_p),form=trim(rstsave),status=stat)    
		write (runit) mpg(1),mpg(2),mpg(3)
		write (runit) p
		close (runit, status='keep')

#if 0
     if(iscalon.eq.1) then

        CALL filename_GEN(filename_scal,'scal',count_routine)
        
        open(unit=runit,file=TRIM(filename_scal),form=trim(rstsave),status=stat)    
        write(runit) nu_error, nu_old
        write(runit) nu_error_array(1:nerr_steps)
#if PARALLEL
        write(runit)phif(0:nx+1,1:my2,1:mz,1:nspmx)
#else
        write(runit)phif
#endif
        write(runit)nlphif
        write(runit) phisurfall(1:nbody,1:nspmx)
        write(runit) phirmean(1:nspmx), fphirmean(1:nspmx)
        close(runit, status='keep')
     endif
#endif


		if (I_AM_NODE_ZERO) then
			filename = TRIM(RUN_NAME)//'_RESTART'
			WRITE(filename1, '(I1)') count_routine
			open(unit=runit,file=TRIM(filename),form="formatted",status="replace")    
			write(runit, *) count_routine, nproc, nprocy, nprocz
			close(runit, status='keep')
		endif


		IF (count_routine.eq.1) count_tmp = 2
		IF (count_routine.eq.2) count_tmp = 1
!!$    PRINT*,'count_tmp = ', myid, count_tmp
		CALL delete_restart_files(count_tmp)

		if(I_AM_NODE_ZERO)then 
			call separator(ounit,62,'-')
			write (ounit,12,ADVANCE='NO')t
			write (ounit,*) TRIM(filename)
			call separator(ounit,62,'-')
		endif
12		format('Saving restart data at time = ',E12.5,' in file with extension ')
	END SUBROUTINE save_unformatted

	SUBROUTINE delete_restart_files(count_routine)
		INTEGER, intent(in):: count_routine
		CHARACTER*100 :: filename, form, filenameLOC
		CHARACTER*8 :: filename1
		INTEGER :: runitno, node, strlen
		LOGICAL:: filexist, isopen
    
		CALL filename_GEN(filename,'u',count_routine)

		INQUIRE(FILE=filename,EXIST=filexist,OPENED=isopen)
		IF(filexist) THEN 
			runitno = getnewunit(minunitno, maxunitno)
			CALL delete_file(filename,trim(rstsave),runitno)

			CALL filename_GEN(filename,'nl',count_routine)
			CALL delete_file(filename,trim(rstsave),runitno)

#if PARALLEL
			CALL filename_GEN(filename,'neighb',count_routine)
			CALL delete_file(filename,trim(rstsave),runitno)
#endif
			CALL filename_GEN(filename,'p',count_routine)
			CALL delete_file(filename,trim(rstsave),runitno)

			IF(ISCALON.EQ.1) then 
				CALL filename_GEN(filename,'scal',count_routine)             
				CALL delete_file(filename,trim(rstsave),runitno)
			end IF
			IF (IMOVE.EQ.1) then
				CALL filename_GEN(filename,'sphere_config',count_routine)             
				IF (I_AM_NODE_ZERO) then
					CALL delete_file(filename,trim(rstsave),runitno)
				ENDIF
			ENDIF
		ENDIF
  end SUBROUTINE delete_restart_files
  
  SUBROUTINE delete_file(filename,formtype,runitno)
    CHARACTER*100, INTENT(in) :: filename
    CHARACTER(LEN=11), INTENT(in):: FORMTYPE
    INTEGER, Intent(in) :: runitno
!    PRINT*,'*******', myid, filename
    OPEN(runitno, FILE=trim(filename), form=formtype, status="replace")
    close(runitno, status='delete')
  end SUBROUTINE delete_file
  

  !------
  ! Read from restart file
  !-------
  
	SUBROUTINE read_restart
		use init_turb, only : forced_turbulence, epsf_forcing, sampling_dissip_time, eddy_time_i
		Implicit none
		character*100 filename
		logical :: filexist

		!INQUIRE(file=restartfile, exist=fileexists)
		!If (fileexists) then
		If (rstread.eq."formatted") then
			! call read_formatted
		elseif (rstread.eq."unformatted") then
			call read_unformatted
		Endif
		!else
		if(I_AM_NODE_ZERO)then
			write(*,12) t
12			format('Read field from restartfile at time = ', F12.5)
		endif

!if (I_AM_NODE_ZERO) then
!write (*,*) "ADSFSDFSDFSDFSDFSDFSDF"
!write (*,*) iturbon, forced_turbulence
!write (*,*) t
!write (*,*) sampling_dissip_time*eddy_time_i
!write (*,*) sampling_dissip_time
!write (*,*) eddy_time_i
!write (*,*) "ADSFSDFSDFSDFSDFSDFSDF"
!endif



		if (iturbon.and.forced_turbulence .and. t>sampling_dissip_time*eddy_time_i) then
			if (I_AM_NODE_ZERO) then
				filename = trim(run_name)//"_linear_forcing.dat"
				inquire (file=trim(filename), exist=filexist)
				if (filexist) then
					open (unit=1, file=trim(filename), status="old", action="read")
					read (1,*) epsf_forcing
					close(1)

					write (*,"(1a,1d15.7)") "DISSIPATION TO CALCULATE FORCING TERM = ", epsf_forcing
				endif
			endif
			BROADCAST_LOGICAL(filexist,1,NODE_ZERO,comm_cart_2d)
			BROADCAST_DOUBLE(epsf_forcing,1,NODE_ZERO,comm_cart_2d)

			if (.not.filexist) then
				if (I_AM_NODE_ZERO) write (*,*) "FORCING TERM FILE NOT FOUND"
				PARALLEL_FINISH()
				stop
			endif
		endif
	END SUBROUTINE read_restart
  
	SUBROUTINE read_unformatted
		use fftw3_interface
		Implicit none
		INTEGER :: count_tmp, i, nerr_steps_tmp, node,strlen, iphs
		!CHARACTER*100 :: filename, filenameLOC

		CHARACTER*100 :: filename_config, filename_neighb, filename_u, filename_p, filename_nl, filename_scal

		CHARACTER*10 :: filename1
		LOGICAL :: filex,isopen

#if PARALLEL
		integer :: error, int_tmp, myid_old, fromy_old, toy_old, &
					&	fromz_old, toz_old, local_ni_old(3), local_no_old(3)
#endif

		if(I_AM_NODE_ZERO)then
			CALL screen_separator(80,'*')
			WRITE(*,*) 'IN READ UNFORMATTED RESTART'
		endif
    
		runit = getnewunit(minunitno,maxunitno)
		if (runit.lt.0) call printerror("newunit","runit")
		count_tmp = count_restart

#if PARALLEL
		CALL filename_GEN(filename_neighb,'neighb',count_tmp)
		open(unit=runit,file=trim(filename_neighb),form="formatted") 
		read(runit,*) myid_old, fromy_old, toy_old, fromz_old, toz_old
		read(runit,*) local_ni_old(1), local_ni_old(2), local_ni_old(3)
		read(runit,*) local_no_old(1), local_no_old(2), local_no_old(3)
		close(runit, status='keep')               

		error = 0
		if (myid_old  /= myid)  error = error+1
		if (fromy_old /= fromy) error = error+1
		if (toy_old   /= toy)   error = error+1
		if (fromz_old /= fromz) error=error+1
		if (toz_old   /= toz)   error= error+1

		if (local_ni_old(1) /= local_ni(1)) error = error+1
		if (local_ni_old(2) /= local_ni(2)) error = error+1
		if (local_ni_old(3) /= local_ni(3)) error = error+1
		if (local_no_old(1) /= local_no(1)) error = error+1
		if (local_no_old(2) /= local_no(2)) error = error+1
		if (local_no_old(3) /= local_no(3)) error = error+1

		if (error>0) then
			write (*,*) "DOMAIN DECOMPOSITION IN DOMAIN : ", myid
		endif


		int_tmp = error
		GLOBAL_INT_SUM(int_tmp, error, 1, comm_cart_2d)
		if (error==0) then
			if (I_AM_NODE_ZERO) write (*,*) "DOMAIN DECOMPOSITION CONSISTENT WITH PREVIOUS TIME..."
		else
			if (I_AM_NODE_ZERO) write (*,*) "DOMAIN DECOMPOSITION NOT CONSISTENT WITH PREVIOUS TIME!!!"

			PARALLEL_FINISH()
			stop
		endif
#endif

		CALL filename_GEN(filename_u,'u',count_tmp)
		open(unit=runit,file=trim(filename_u),form=trim(rstsave))
    
		read(runit) tstart, dt
		read(runit)umean(1), umean(2), umean(3)
		read(runit)usmean_des(1), usmean_des(2), usmean_des(3)
		read(runit)ufmean_des(1), ufmean_des(2), ufmean_des(3)
		read(runit)cf
    
		t = tstart

		cforig = cf
		if (I_AM_NODE_ZERO) Write(*,*)'cf from restart = ', cf
		read(runit) nerr_steps_tmp

		read(runit)ferror, fold

		if (SIZE(ferror_array,1).ne.nerr_steps) then
			DEALLOCATE(ferror_array)
			ALLOCATE(ferror_array(nerr_steps))
		endif
    
		if(nerr_steps_tmp.gt.nerr_steps) then
			read(runit) ferror_array(1:nerr_steps)
		else
			ferror_array = 1.d0
			read(runit) ferror_array(1:nerr_steps_tmp)
		end if

		do iphs = 1, nphases
			read(runit)phase_array(iphs)%ferror, phase_array(iphs)%fold
			!^^^^^^^ 06-24-2010 Mohammad: RESTARTING A FIXED BED FOR MOVING CASE ^^^^
			if (imove==1) then
				if (moving_from_fixed_bed) then
					phase_array(iphs)%grant_error = one
					phase_array(iphs)%grant_old = gran_temp
				else
					read(runit)phase_array(iphs)%grant_error, phase_array(iphs)%grant_old
				endif
			endif
			!------------------------------------------------------------------------
		enddo

		do iphs = 1, nphases
			if(SIZE(phase_array(iphs)%ferror_array,1).ne.nerr_steps)then
				if(ASSOCIATED(phase_array(iphs)%ferror_array)) DEALLOCATE(phase_array(iphs)%ferror_array)
				ALLOCATE(phase_array(iphs)%ferror_array(nerr_steps))
			end if
			if(imove.eq.1)then
				if(SIZE(phase_array(iphs)%grant_array,1).ne.nerr_steps)then
					DEALLOCATE(phase_array(iphs)%grant_array)
					ALLOCATE(phase_array(iphs)%grant_array(nerr_steps))
				end if
			end if
		end do

		do iphs = 1, nphases
			if(nerr_steps_tmp.gt.nerr_steps) then
				read(runit) phase_array(iphs)%ferror_array(1:nerr_steps)
			else
				phase_array(iphs)%ferror_array = 1.d0
				read(runit) phase_array(iphs)%ferror_array(1:nerr_steps_tmp)
			end if
			if(imove.eq.1)then
				!^^^^^^^ 06-24-2010 Mohammad: RESTARTING A FIXED BED FOR MOVING CASE ^^^^          
				if (moving_from_fixed_bed) then
					phase_array(iphs)%grant_array(1:nerr_steps) = one
				else
					if(nerr_steps_tmp.gt.nerr_steps) then
						read(runit) phase_array(iphs)%grant_array(1:nerr_steps)
					else
						phase_array(iphs)%grant_array = 1.d0
						read(runit) phase_array(iphs)%grant_array(1:nerr_steps_tmp)
					end if
				endif
				!--------------------------------------------------------------------------
			end if
		end do

		read (runit) u
		close (runit, status='keep')
    
		if (I_AM_NODE_ZERO) then
			WRITE(*,'(A40,g12.5)')'READING THE RESTART FILE WRITTEN AT ', tstart
			WRITE(*,'(A40,g12.5)')'DT READ FROM RESTART FILES = ', dt
		endif

		ferror_hist = SUM(ferror_array(1:nerr_steps))/nerr_steps
		if(I_AM_NODE_ZERO)then
			WRITE(*,'(A25,2(2x,g12.5))')'FOLD, FERROR = ', FOLD, FERROR
			WRITE(*,'(A25,g12.5)') 'FERROR_HIST = ', ferror_hist
		end if
    
		do iphs = 1, nphases
			phase_array(iphs)%ferror_hist = SUM(phase_array(iphs)%ferror_array(1:nerr_steps))/nerr_steps
			if(I_AM_NODE_ZERO)then
				WRITE(*,'(A10,I4,A25,2(2x,g12.5))')'PHASE = ', iphs, 'FOLD, FERROR = ', phase_array(iphs)%fold, phase_array(iphs)%ferror
				WRITE(*,'(A10,I4,A25,g12.5)') 'PHASE = ', iphs, 'FERROR_HIST = ', phase_array(iphs)%ferror_hist
			endif
		enddo
    
		CALL filename_GEN(filename_nl,'nl',count_tmp)
		open(unit=runit,file=trim(filename_nl),form=trim(rstsave))
		read (runit) nl
		close(runit, status='keep')

		CALL filename_GEN(filename_p,'p',count_tmp)
		open(unit=runit,file=trim(filename_p),form=trim(rstsave))
		read(runit)mpg(1),mpg(2),mpg(3)
		read(runit)p
		close(runit, status='keep')
    
#if 0
    if(iscalon.eq.1.and.iscal_restart.eq.1) then
       CALL filename_GEN(filename_scal,'scal',count_tmp)
       open(unit=runit,file=trim(filename_scal),form=trim(rstsave))
       
       read(runit) nu_error, nu_old
       
       if(nerr_steps_tmp.gt.nerr_steps) then
          read(runit) nu_error_array(1:nerr_steps)
       else
          nu_error_array = 1.d0
          read(runit) nu_error_array(1:nerr_steps_tmp)
       end if
#if PARALLEL
       read(runit)phif(0:nx+1,:,:,1:nspmx)
#else
       read(runit)phif
#endif
       read(runit) nlphif(1:nx,:,:,1:nspmx)
       read(runit) phisurfall(1:nbody,1:nspmx)
       read(runit) phirmean(1:nspmx),  fphirmean(1:nspmx)
       close(runit, status = "keep")
       
       if(I_AM_NODE_ZERO)then
          WRITE(*,*) 'IN SCALAR READ UNFORMATTED RESTART'
          
          DO i=1,nbody
             !WRITE(*,'(A,g12.5,2x,i2)')'phisurf of bodies =', phisurfall(i,1),i
             WRITE(ounit,'(A,g12.5,2x,i2)')'phisurf of bodies =', phisurfall(i,1),i
          END DO
          WRITE(*,'(A30,2(2x,g12.5))')'phirmean and fphirmean = ', phirmean, fphirmean 
          
          WRITE(*,'(A25,2(2x,g12.5))')'nu_old, nu_error = ', nu_old, nu_error
       end if
       nu_error_hist = SUM(nu_error_array(1:nerr_steps))/nerr_steps
       if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.5)')'NU_ERROR_HIST = ', NU_ERROR_HIST
    endif
#endif
		if(I_AM_NODE_ZERO)CALL screen_separator(80,'*')
    
	END SUBROUTINE read_unformatted
end MODULE restart_funcs
