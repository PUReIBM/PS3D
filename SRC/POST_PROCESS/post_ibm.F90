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

PROGRAM post_ibm
#include "../FLO/ibm.h"
  !
  ! AUTHOR: RAHUL GARG 
  ! POST PROCESSOR CODE FOR DTIBM

  USE precision  
  USE restart_funcs
  USE constants  
  USE post_global_data
  USE initialize_flo 
  USE nlcalc
  USE steptheflow
  USE initialize  
  USE writeoutput 
  USE machine 
  USE scalar_post
  USE spectra_3d 
  USE stepthescalar
  USE nlarrays, ONLY : dvf=>uf1, dvr=>ur1, dpf=>uf2, dpr=>ur2
  USE langevin_model
  USE lang_source_diss
  USE useful_functions
  USE scalar_data 
  USE outputscalar 
  USE post_static
  USE post_freely_evolving
  IMPLICIT NONE

  real(prcn) :: cpu0, cpu1, cputime_used, cputime_hrs, cputime_min

  CHARACTER*80 :: post_file, conv_file, filename, blowup_file
  LOGICAL :: filexist,isopen
  CHARACTER::         MIS*80
  INTEGER :: countmis,idim,jdim,iphs, tunit, i, isp, imis, PP_MIS,&
       & misavgdragunit, imeasvol, numvollineunit
  REAL(prcn) :: avg_nusselt, avg_nusselt_para,&
       & avg_nusselt_perp, err_bar_nu, err_bar_nu_para,&
       & err_bar_nu_perp,  avg_solid_area_frac(1:mxmax), err_bar_solid_area_frac(1:mxmax)
  REAL(prcn), DIMENSION(mxmax,3):: avg_scal_frac_diff, err_bar_scal_diff,avg_scal_frac_nl, err_bar_scal_nl,avg_scal_frac_nl_diff, err_bar_scal_nl_diff,avg_scal_frac_fphi, err_bar_scal_fphi

  CHARACTER*80 :: table_name
  LOGICAL :: post_process, blowup_filexist, specified_runname

  PARALLEL_START()
  GET_NPROCS(comm_group,nproc)

  GET_PROCESSOR_RANK(comm_group,myid)

  specified_runname = .FALSE.
  RUN_NAME = ""
  if(I_AM_NODE_ZERO)then
     CALL GETARG(1 , POST_RUNNAME)
     CALL GETARG(2, MIS)
     READ(MIS, FMT='(I5)') nmis
     if(nmis.eq.1)then
        CALL GETARG(3, RUN_NAME)
     end if
     !CALL GETARG(4, PER_CONF)
     if(TRIM(RUN_NAME).eq."")then
        SPECIFIED_RUNNAME = .FALSE.
     else
        SPECIFIED_RUNNAME = .TRUE.
     end if
  end if
  BROADCAST_INT(nmis, 1, NODE_ZERO,comm_group)
  BROADCAST_LOGICAL(SPECIFIED_RUNNAME, 1, NODE_ZERO,comm_group)

  IF(nmis.gt.nmismax) then
     if(I_AM_NODE_ZERO)WRITE(*,'(A,\,A)') "WARNING: nmis GT nmismax used in post_global_data"," STOPPING"
     PARALLEL_FINISH()
     STOP
  end IF

  !if(TRIM(per_conf).eq."")
  per_conf='95'

  PER_CONF = TRIM(PER_CONF)//'%'
  PRINT*,'PERCENTAGE CONFIDENCE = ', PER_CONF
  CALL calculate_constants
  CALL GET_RUN_ID 
  WRITE(*, *) 'Post processing started on date:', ID_MONTH, ID_DAY, ID_YEAR, ' AT TIME:', ID_HOUR, ID_MINUTE , ID_SECOND
  WRITE(*,*) 'RUNNING ON MACHINE ', ID_NODE
  CALL CPU_TIME (CPU0) 

  FROM_POST = .true. 
  POST_NO_FLOW_MEM_ALLOC = .FALSE.
  countmis = 0

  tke_fluid = zero
  tke_fluid_var = zero
  MIS_CONV(:) = 0

  !NULLIFY(mis_data)
  !NULLIFY(current_mis)
  !ALLOCATE(mis_data)
  !current_mis => mis_data
!!$  !NULLIFY(mis_data%next)

  ENSAVG%force_max(1:ndim) = -SMALL_NUMBER
  ENSAVG%force_min(1:ndim) = LARGE_NUMBER

  ENSAVG%vel_max(1:ndim) = -SMALL_NUMBER
  ENSAVG%vel_min(1:ndim) = LARGE_NUMBER

  do imis = 1,nmis
     Write(*,*) 'Specified RUNNAME = ', specified_runname
     if(.not.specified_runname)then
        if(imis.lt.10)then
           WRITE(RUN_NAME,fmt="('MIS',i1)")imis
        else
           WRITE(RUN_NAME,fmt="('MIS',i2)")imis
        end if
     end if
     post_file = TRIM(RUN_NAME)//"_POST_PROCESSED"
     conv_file = TRIM(RUN_NAME)//"_CONVERGED"
     blowup_file = TRIM(RUN_NAME)//"_BLOWUP.dat"
     INQUIRE(FILE=conv_file,EXIST=filexist,OPENED=isopen)
     !INQUIRE(FILE=TRIM(RUN_NAME)//'_floparam.in',EXIST=filexist,OPENED=isopen)
     post_process = filexist
     post_process = .TRUE.
     INQUIRE(FILE=blowup_file,EXIST=blowup_filexist,OPENED=isopen)

     MISTRUE: if(post_process.and.(.not.blowup_filexist))then
        if(ALLOCATED(phase_array))CALL DEALLOC_PHASE_RELATED
        if(ALLOCATED(xc))  CALL DEALLOC_BND_RELATED
        !Set mis_conv to true 
        MIS_CONV(imis) = 1
        countmis = countmis +1
        CALL screen_separator(100, '-')

        CALL screen_separator(100, 'C')
        if(imove.eq.1)then 
           WRITE(*,'(A)')'THIS IS A MOVING PARTICLE CASE. POST PROCESSING ANYWAYS'
        else
           WRITE(*,'(A,A)')'CONGRATULATIONS: CONVERGED FILE EXISTS FOR RUN NAME = ', RUN_NAME
           WRITE(*,'(A,i4)') 'NUMBER OF CONVERGED MIS SO FAR =', SUM(MIS_CONV)
        end if

        CALL screen_separator(100, '-')

        CALL screen_separator(100, 'C')

        CALL init_params !set the parameters

        irestart = 1
        iscal_restart = irestart 
        flow_converged = .true.
!!$        nerr_steps = 50
!!$        tol_ferror = 1.0E-05 
        CALL initflo
        WRITE(*,*) 'DONE WITH INITIALIZATION. CALLING NON LINEAR'
        FIRST_PASS = .FALSE.
        
        WRITE(*,*) 'WILL NOW POST PROCESS THIS CASE'

        !CALL post_mean_drag_poly

        !IF(RET.GT.ZERO)then
        !CALL populate_current_realization(imis)
        !CALL calc_source_diss
        !end IF

        !CALL calc_spec_3d

        if(imove.eq.1)then
           if(.not.move_particles) then
              WRITE(*,'(A)')'PARTILCES HAVE NOT STARTED MOVING YET. EXITING.'
              PARALLEL_FINISH()
              STOP
           end if
           ! MOVING PARTICLES CASE. CALL THE RELEVANT POST PROCESSING FUNCTION
           WRITE(*,'(A)')'POST PROCESSING MOVING PARTICLE CASE.'
           CALL post_moving_particles(imis)
        else
           CALL post_static_particles(imis)
        end if

        if(iscalon.eq.1)then 

           IF(mx.gt.mxmax) then
              WRITE(*,'(A,\,A)') "WARNING: mx GT mxmax used in post_global_data"," MIGHT GET SEG SEV ERRORS LATER"
           end IF
           WRITE(*,*) 'nspmx at file reading = ', nspmx
           CALL output_scal 
           CALL scalar_pp(imis)
           !           WRITE(*,*) 'nspmx after scalar pp = ', nspmx
        end if


        OPEN(unit = 1000, file=post_file, status="replace")
        close(1000, status="keep")

     ELSE 

        CALL screen_separator(100, '-')

        CALL screen_separator(100, 'S')
        WRITE(*,'(A,A)')'CONVERGED FILE DOES NOT EXIST FOR RUN NAME = ', RUN_NAME

        CALL screen_separator(100, '-')

        CALL screen_separator(100, 'S')

     end if MISTRUE

  end do


  if(imove.ne.1)then

     misavgdragunit = getnewunit(minunitno,maxunitno)
     OPEN(misavgdragunit,FILE=TRIM(POST_RUNNAME)//'_mis_avg_drag.dat', form='formatted')
    
     do iphs = 1, nphases
        CALL MIS_AVG_DATA(MIS_NORM_DRAG_SPEC(:,iphs),&
             & AVG_NORM_DRAG_SPEC(iphs), NORM_DRAG_SPEC_ERRBAR(iphs),&
             & confin, mis_conv(:), check_for_outlier = .false.)

        CALL MIS_AVG_DATA(MIS_NORM_DRAG_CHEM_SPEC(:,iphs),&
             & AVG_NORM_DRAG_CHEM_SPEC(iphs), NORM_DRAG_CHEM_SPEC_ERRBAR(iphs),&
             & confin, mis_conv(:), check_for_outlier = .false.)
        Write(misavgdragunit,'(10(2x,g20.12))') phase_array(iphs)%dia/char_length, avg_norm_drag_spec(iphs), avg_norm_drag_chem_spec(iphs), NORM_DRAG_SPEC_ERRBAR(iphs), NORM_DRAG_CHEM_SPEC_ERRBAR(iphs)
     end do
     CLOSE(misavgdragunit,status='keep') 


     numvollineunit = getnewunit(minunitno, maxunitno)
     OPEN(numvollineunit, FILE=TRIM(POST_RUNNAME)//'_NUMVOLLINEPLOT.dat',&
          & form='formatted')
     do imeasvol = 1, nmeasvols
        CALL MIS_AVG_DATA(NUMBER_IN_MIS(:,imeasvol),&
             & NUMBER_IN_MEASVOL(imeasvol), NUMBER_STDDEV(imeasvol),&
             & confin, mis_conv(:), check_for_outlier = .false.)

        CALL MIS_AVG_DATA(NUMBERSQ_IN_MIS(:,imeasvol),&
             & NUMBERSQ_IN_MEASVOL(imeasvol), NUMBERSQ_STDDEV(imeasvol),&
             & confin, mis_conv(:), check_for_outlier = .false.)

        CALL MIS_AVG_DATA(VOLUME_IN_MIS(:,imeasvol),&
             & VOLUME_IN_MEASVOL(imeasvol), VOLUME_STDDEV(imeasvol),&
             & confin, mis_conv(:), check_for_outlier = .false.)

        CALL MIS_AVG_DATA(VOLUMESQ_IN_MIS(:,imeasvol),&
             & VOLUMESQ_IN_MEASVOL(imeasvol), VOLUMESQ_STDDEV(imeasvol),&
             & confin, mis_conv(:), check_for_outlier = .false.)

        Write(numvollineunit,'(20(2x,g17.5))')LMEAS_VOLS(imeasvol)/dia_phys,&
             & NUMBER_IN_MEASVOL(imeasvol)&
             &,NUMBER_STDDEV(imeasvol), &
             & NUMBERSQ_IN_MEASVOL(imeasvol),&
             & NUMBERSQ_STDDEV(imeasvol),&
             & VOLUME_IN_MEASVOL(imeasvol),&
             & VOLUME_STDDEV(imeasvol),&
             & VOLUMESQ_IN_MEASVOL(imeasvol),&
             & VOLUMESQ_STDDEV(imeasvol)
     end do
     CLOSE(numvollineunit, status='keep')
  end if

  
  !IF(RET.GT.ZERO)
  Write(*,*)'WRITING ENSEMBLE AVERAGE QUANTITIES.'
  !CALL write_ensemble_average_data

end PROGRAM post_ibm

  
