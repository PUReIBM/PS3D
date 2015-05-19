MODULE post_freely_evolving
#include "../FLO/ibm.h"
  !
  ! AUTHOR: Sudheer Tenneti
  ! POST PROCESSOR CODE FOR FREELY EVOLVING SUSPENSIONS
  USE precision  
  USE restart_funcs
  USE general_funcs
  USE postproc_funcs
  USE constants  
  USE randomno
  USE useful_functions
  USE dependent_functions
  USE post_global_data

  IMPLICIT NONE

  PRIVATE
  REAL(prcn),ALLOCATABLE, DIMENSION(:) :: granular_temperature,&
       & char_force, meanforcemod, phase_mass, acclnstddev_mod, vivjzioft,&
       & vivjetaoft, aivjzioft, aivjetaoft, source_grantemp&
       &,source_pres_grantemp, source_visc_grantemp, dissip_grantemp,&
       & dissip_pres_grantemp, dissip_visc_grantemp, steady_granular_temp,&
       & steady_source_grantemp, steady_dissip_grantemp, nondim_grantemp,&
       & nondim_source_grantemp, nondim_dissip_grantemp
  REAL(prcn), ALLOCATABLE, DIMENSION(:) :: RAD_PHYS
  REAL(prcn),ALLOCATABLE, DIMENSION(:,:) :: fluctv,fluct_force,&
       & meanvel,varvel,meanforce, meanaccln, fluct_accln,&
       & acclnstddev, fluct_pres_accln, fluct_visc_accln, XC_PHYS
  
  REAL(prcn),ALLOCATABLE, DIMENSION(:,:,:) :: AIVJ, VIVJ

  REAL(prcn) :: currtime, eigenvalue(ndim), eigen_vectors(ndim,ndim),&
       & numdens

  
  TYPE GRID
     REAL(prcn) :: dx, dy, dz
     INTEGER :: cx, cy, cz ! Number of grid nodes
     REAL(prcn), DIMENSION(:), POINTER :: XE, YN, ZT
     REAL(prcn), DIMENSION(:,:,:), POINTER :: meannum, meannumsq
     REAL(prcn), DIMENSION(:,:,:), POINTER :: meanvol, meanvolsq
  END TYPE GRID
  !TYPE(GRID), DIMENSION(:), ALLOCATABLE, TARGET :: coarse_grids
  TYPE(GRID) :: coarse_grid
  LOGICAL ::  grant_converged

  PUBLIC ::  post_moving_particles

CONTAINS

  SUBROUTINE post_moving_particles(imis)
    USE constants  
    IMPLICIT NONE
    INTEGER, Intent(in) :: imis
    INTEGER :: iostatus, idim, iphs, npart, pstart, pend, m, jdim,&
         & timestep_count
    CHARACTER(LEN=80) :: FILENAME
    REAL(prcn) :: symmaivj(nphases,ndim,ndim), maxeigen, dummy1,&
         & dummy2, mphase, diaphase, source_char(nphases),&
         & dissip_char(nphases), source_dissip_char(nphases)

    if(ALLOCATED(granular_temperature))then
       CALL DEALLOCATE_MOVING_PART_RELATED
    end if

    CALL ALLOCATE_MOVING_PART_RELATED
    
    !CALL compute_geometric_statistics
    
    
    FILENAME = TRIM(RUN_NAME)//'_part_info'//'.rst'
    OPEN(1050, FILE = FILENAME, STATUS='old', form='unformatted')
    READ(1050)movestart_time
    CLOSE(1050, status = 'keep')
    
    OPEN(1050, FILE = FILENAME, STATUS='old', form='unformatted')
    WRITE(*,'(A)')'NOW OPENING THE RELEVANT FILES FOR WRITING .'
    CALL open_relevant_files
    
    CALL compute_measvol_stats
    STOP

#if 0    
    OPEN(1051, FILE = TRIM(RUN_NAME)//'_sphr_center_out'//'.dat',&
         & STATUS='old', form='formatted')
    READ(1051, *)
    do m = 1, nbody
       READ(1051,*)xc(m,1), xc(m,2), xc(m,3), dummy1, dummy2, dummy1, dummy2
    end do
    CLOSE(1051,status='keep')
    CALL calculate_gofr_homog
    STOP
#endif

    do iphs = 1, nphases
       char_force(iphs) = 3.d0*pi*vis*(ucharmod)*phase_array(iphs)%dia
       phase_mass(iphs) = RHOS*pi*(phase_array(iphs)%dia)**3.d0/6.d0
    end do
    WRITE(*,'(A)')'NOW POST PROCESSING .'
    grant_converged = .FALSE.
    steady_granular_temp = zero
    steady_source_grantemp = zero
    steady_dissip_grantemp = zero
    timestep_count = 0
    
    Write(*,*)'READING THE FILE...'
125 continue    
    
    READ(1050,END=225)post_proc_time
    READ(1050,END=225)xc(1:nbody, 1:3)
    READ(1050,END=225)velbdy(1:nbody, 1:3)
    READ(1050,END=225)force(1:nbody, 1:3)
    READ(1050,END=225)pres(1:nbody, 1:3)
    READ(1050,END=225)visc(1:nbody, 1:3)
    READ(1050,END=225) frame_vel(1:3)
    READ(1050,END=225) frame_accln(1:3)
    READ(1050,END=225) ufmean(1:3)
    
    currtime = post_proc_time - movestart_time
    

    ! Compute the principle axes of moment of inertia for the whole
    ! suspension
    !CALL compute_moment_of_inertia_principal_coordinates

    !CALL calculate_gofr_homog
    ! Compute the mean and fluctuating quantites

    
    CALL compute_mean_and_fluct_quantites

    ! Compute the granular temperature and start storing the time trace in a file
    
    CALL compute_granular_temperature(.TRUE.)

!!$    if(.not.grant_converged)then
!!$       grant_converged = .TRUE.
!!$       do iphs = 1, nphases
!!$          grant_converged = grant_converged.and.(phase_array(iphs)%gran_error_hist.lt.tol_gran_error)
!!$       end do
!!$    end if
    
    
    do iphs = 1, nphases
       npart = phase_array(iphs)%npart
       pstart = phase_array(iphs)%pstart
       pend = phase_array(iphs)%pend

       CALL compute_grantemp_source_dissipation(iphs)       

       CALL compute_aivj(aivj(iphs,1:ndim,1:ndim)&
            &,fluct_accln(pstart:pend,1:ndim),fluctv(pstart:pend&
            &,1:ndim),npart)

       CALL compute_aivj(vivj(iphs,1:ndim,1:ndim)&
            &,fluctv(pstart:pend,1:ndim),fluctv(pstart:pend&
            &,1:ndim),npart)
       do idim = 1, ndim
          do jdim = 1, ndim
             symmaivj(iphs,idim,jdim) = aivj(iphs,idim,jdim) +&
                  & aivj(iphs,jdim,idim)
          end do
       end do
       CALL calc_anisotropy(symmaivj(iphs,1:ndim,1:ndim), aivjzioft(iphs),&
            & aivjetaoft(iphs))
       CALL calc_anisotropy(vivj(iphs,1:ndim,1:ndim), vivjzioft(iphs),&
            & vivjetaoft(iphs))
    end do

    WRITE(aivjunit,'(500(2x, e20.12))') currtime, currtime/t_conv&
         &, currtime/t_vis, (((aivj(iphs,idim,jdim)/(acclnstddev(iphs,idim)&
         &*DSQRT(varvel(iphs,idim))),jdim=1,ndim),idim=1,ndim)&
         &,iphs=1,nphases)
    !WRITE(anisunit,*)'ZONE T = " ', currtime/t_conv, ' "'
    WRITE(anisunit, '(10(2x, g17.5))')currtime/t_conv, (vivjzioft(iphs), vivjetaoft(iphs),&
         & aivjzioft(iphs), aivjetaoft(iphs), iphs = 1, nphases)

    do iphs = 1, nphases
       diaphase = phase_array(iphs)%dia
       mphase = pi*diaphase**3.d0/6.d0
       source_char(iphs) = (3.d0*pi*vis*ucharmod*diaphase/mphase) *&
            & DSQRT(granular_temperature(iphs)+SMALL_NUMBER)
       dissip_char(iphs) = (3.d0*pi*vis*diaphase/mphase)&
            &*(granular_temperature(iphs)+SMALL_NUMBER)
       source_dissip_char(iphs) = 3.d0*pi*vis*diaphase*ucharmod**2.d0/mphase
    end do
    
    WRITE(source_dissip_unit, '(20(2x&
         &,g17.5))')currtime/t_conv, (granular_temperature(iphs),&
         & source_grantemp(iphs), dissip_grantemp(iphs),iphs = 1&
         &,nphases), (source_grantemp(iphs)/source_char(iphs),&
         & dissip_grantemp(iphs)/dissip_char(iphs),&
         & source_grantemp(iphs)/source_dissip_char(iphs)&
         &,dissip_grantemp(iphs)/source_dissip_char(iphs)&
         &,(source_grantemp(iphs)-dissip_grantemp(iphs))&
         &/source_dissip_char(iphs), iphs = 1, nphases)

    WRITE(source_dissip_pres_visc_unit, '(9(2x&
         &,g17.5))')currtime/t_conv, (granular_temperature(iphs),&
         & source_grantemp(iphs), dissip_grantemp(iphs),&
         & source_pres_grantemp(iphs), dissip_pres_grantemp(iphs),&
         & source_visc_grantemp(iphs), dissip_visc_grantemp(iphs),&
         & iphs = 1,nphases)



!!$    if(grant_converged)then
!!$       do iphs = 1, nphases
!!$          steady_granular_temp(iphs) = steady_granular_temp(iphs) + granular_temperature(iphs)
!!$          steady_source_grantemp(iphs) = steady_source_grantemp(iphs) + source_grantemp(iphs)
!!$          steady_dissip_grantemp(iphs) = steady_dissip_grantemp(iphs) + dissip_grantemp(iphs)
!!$       end do
!!$       timestep_count = timestep_count + 1
!!$    end if

    goto 125

225 continue

    Write(*,*)'REACHED END OF FILE. CLOSING THE FILE'
    CALL close_relevant_files
    CLOSE(1050, status = 'keep')    
    
!!$    steady_granular_temp(1:nphases) = steady_granular_temp(1:nphases)/real(timestep_count,prcn)
!!$    steady_source_grantemp(1:nphases) = steady_source_grantemp(1:nphases)/real(timestep_count,prcn)
!!$    steady_dissip_grantemp(1:nphases) = steady_dissip_grantemp(1:nphases)/real(timestep_count,prcn)
    !Write(*,*)'Steady granular temperature = ', steady_granular_temp(1:nphases)
    !CALL compute_steady_state_statistics
    
  END SUBROUTINE post_moving_particles



  SUBROUTINE compute_mean_and_fluct_quantites
    IMPLICIT NONE
    INTEGER :: idim, iphs, npart, pstart, pend, m
    REAL(prcn) :: u(nbody),wt(nbody),adev,var,skew,curt,sdev
    REAL(prcn) :: meanpresaccln(nphases,ndim), meanviscaccln(nphases,ndim)

    do iphs = 1, nphases
       npart = phase_array(iphs)%npart
       pstart = phase_array(iphs)%pstart
       pend = phase_array(iphs)%pend

       ! Mean and fluctuating velocity
       do idim=1, ndim
          u(1:npart) = velbdy(pstart:pend,idim)
          wt(1:npart) = 1/real(npart,prcn)
          CALL moment1(4, npart, npart, wt(1:npart), u(1:npart),&
               & meanvel(iphs,idim),adev,sdev,varvel(iphs,idim),skew,curt)

          do m = pstart, pend
             fluctv(m,idim) = (velbdy(m,idim)-meanvel(iphs,idim))
          end do
       end do

       ! Mean and fluctuating force
       do idim=1, ndim
          u(1:npart) = force(pstart:pend,idim)
          wt(1:npart) = 1/real(npart,prcn)
          CALL moment1(4, npart, npart, wt(1:npart), u(1:npart),&
               & meanforce(iphs,idim),adev,sdev,var,skew,curt)

          do m = pstart, pend
             fluct_force(m,idim) = (force(m,idim)-meanforce(iphs,idim))
          end do
       end do

       meanforcemod(iphs) = DOT_PRODUCT(meanforce(iphs,1:ndim),meanforce(iphs,1:ndim))
       meanforcemod(iphs) = DSQRT(meanforcemod(iphs))


       ! Mean and fluctuating acceleration
       do idim=1, ndim
          u(1:npart) = force(pstart:pend,idim)/phase_mass(iphs)
          wt(1:npart) = 1/real(npart,prcn)
          CALL moment1(4, npart, npart, wt(1:npart), u(1:npart),&
               & meanaccln(iphs,idim),adev,acclnstddev(iphs,idim),var&
               &,skew,curt)

          do m = pstart, pend
             fluct_accln(m,idim) = (force(m,idim)/phase_mass(iphs)-meanaccln(iphs,idim))
          end do
       end do
       acclnstddev_mod(iphs) = DOT_PRODUCT(acclnstddev(iphs,1:ndim)&
            &,acclnstddev(iphs,1:ndim))
       acclnstddev_mod(iphs) = DSQRT(acclnstddev_mod(iphs))

       ! Mean and fluctuating acceleration from fluctuating pressure
       do idim=1, ndim
          u(1:npart) = pres(pstart:pend,idim)/phase_mass(iphs)
          wt(1:npart) = 1/real(npart,prcn)
          CALL moment1(4, npart, npart, wt(1:npart), u(1:npart),&
               & meanpresaccln(iphs,idim),adev,sdev,var&
               &,skew,curt)

          do m = pstart, pend
             fluct_pres_accln(m,idim) = (pres(m,idim)/phase_mass(iphs)-meanpresaccln(iphs,idim))
          end do
       end do

       ! Mean and fluctuating acceleration from viscous contribution
       do idim=1, ndim
          u(1:npart) = visc(pstart:pend,idim)/phase_mass(iphs)
          wt(1:npart) = 1/real(npart,prcn)
          CALL moment1(4, npart, npart, wt(1:npart), u(1:npart),&
               & meanviscaccln(iphs,idim),adev,sdev,var&
               &,skew,curt)

          do m = pstart, pend
             fluct_visc_accln(m,idim) = (visc(m,idim)/phase_mass(iphs)-meanviscaccln(iphs,idim))
          end do
       end do

!!$       Write(*,*) 'PRES + VISC = ', fluct_pres_accln(1,1)&
!!$            &+fluct_visc_accln(1,1), 'FLUCT ACCLN = ', fluct_accln(1,1)
    end do

  END SUBROUTINE compute_mean_and_fluct_quantites

  SUBROUTINE compute_granular_temperature(IFWRITE)
    IMPLICIT NONE
    LOGICAL, Intent(in) :: IFWRITE
    INTEGER :: iphs, npart, pstart, pend, idim
    REAL(prcn) :: granT(nphases), temp_grant_array(nerr_steps),&
         & ret_curr(nphases), T_par(nphases), T_per(nphases)
    LOGICAL, SAVE :: FIRST_TIME_HERE = .TRUE.

    do iphs = 1, nphases
       npart = phase_array(iphs)%npart
       pstart = phase_array(iphs)%pstart
       pend = phase_array(iphs)%pend

       granT(iphs) = zero
       granT(iphs) = varvel(iphs,1) + varvel(iphs,2) + varvel(iphs,3)
       grant(iphs) = one/three * grant(iphs)
       granular_temperature(iphs) = grant(iphs)

       if(IFWRITE)then
          if(first_time_here)then
             phase_array(iphs)%grant_old = grant(iphs)
             first_time_here = .FALSE.
          end if

          if(grant(iphs).gt.zero)then
             phase_array(iphs)%grant_error = ABS(phase_array(iphs)%grant_old - grant(iphs))/grant(iphs)
          else
             phase_array(iphs)%grant_error = one
          end if
          phase_array(iphs)%grant_old = grant(iphs)

          temp_grant_array(1:nerr_steps) = phase_array(iphs)%grant_array(1:nerr_steps)
          phase_array(iphs)%grant_array(2:nerr_steps) = temp_grant_array(1:nerr_steps-1)
          phase_array(iphs)%grant_array(1) = phase_array(iphs)%grant_error

          phase_array(iphs)%gran_error_hist = SUM(phase_array(iphs)&
               &%grant_array(1:nerr_steps))/nerr_steps

          Ret_curr(iphs) = phase_array(iphs)%dia * DSQRT(grant(iphs))/vis
!!$       Write(*,'(A25,2x,I,2x,A2,2x,g12.5)')'RE_T OF PHASE ', iphs,
          !! ' = ', Ret_curr(iphs)
          T_par(iphs) = varvel(iphs,1)
          T_per(iphs) = half*(varvel(iphs,2)+varvel(iphs,3))
         
       end if
    end do

    if(IFWRITE)then
       WRITE(grantempunit,'(500(2x, e20.12))') currtime, currtime/t_conv&
            &, currtime/t_vis, granular_temperature(1:nphases)&
            &/(ucharmod**2.d0),(Ret_curr(iphs), iphs = 1, nphases)&
            &,(T_par(iphs)/ucharmod**2.d0, iphs = 1, nphases),&
            & (T_per(iphs)/(ucharmod**2.d0), iphs = 1, nphases),&
            & (T_par(iphs)/T_per(iphs),iphs = 1, nphases),&
            & (phase_array(iphs)%gran_error_hist, iphs = 1, nphases),&
            & ufmean(1:ndim), meanvel(1,1:ndim), frame_vel(1:ndim),&
            & meanforcemod(1)/char_force(1)
    end if

  END SUBROUTINE compute_granular_temperature

  SUBROUTINE compute_grantemp_source_dissipation(iphs)
    IMPLICIT NONE
    INTEGER, Intent(in) :: iphs
    INTEGER :: idim, npart, pstart, pend, m
    REAL(prcn) :: fluctvmag2
    REAL(prcn) :: zetan, zetanpres, zetanvisc
    REAL(prcn) :: abszetan, abszetanpres, abszetanvisc

    npart = phase_array(iphs)%npart
    pstart = phase_array(iphs)%pstart
    pend = phase_array(iphs)%pend

    source_grantemp(iphs) = zero
    dissip_grantemp(iphs) = zero
    source_pres_grantemp(iphs) = zero
    dissip_pres_grantemp(iphs) = zero
    source_visc_grantemp(iphs) = zero
    dissip_visc_grantemp(iphs) = zero

    do m = pstart, pend
       fluctvmag2 = DOT_PRODUCT(fluctv(m,1:ndim), fluctv(m, 1:ndim))

       zetan = -DOT_PRODUCT(fluct_accln(m,1:ndim), fluctv(m, 1:ndim))
       zetan = zetan/fluctvmag2
       abszetan = ABS(zetan)
       !if(zetan.gt.zero)Write(*,*)'zetan is positive for m = ', m
       zetanpres = -DOT_PRODUCT(fluct_pres_accln(m,1:ndim), fluctv(m, 1:ndim))
       zetanpres = zetanpres/fluctvmag2
       abszetanpres = ABS(zetanpres)

       zetanvisc = -DOT_PRODUCT(fluct_visc_accln(m,1:ndim), fluctv(m, 1:ndim))
       zetanvisc = zetanvisc/fluctvmag2
       abszetanvisc = ABS(zetanvisc)
       !if(zetanvisc.gt.zero)then
       !Write(*,*)'zetan is positive for m = ', m
       !Write(*,*)'zetanvisc = ', zetanvisc, 'zetanpres = ',&
       !        & zetanpres
       !   Write(*,*)'zetan = ', zetan, 'zetanpres+visc = ',&
       !        & zetanpres+zetanvisc
       !end if

       
!!$       if(m.eq.1)then
!!$          Write(*,*)'zetan = ', zetan, 'zetan pres , zetan visc = ',&
!!$               & zetanpres, zetanvisc
!!$       end if
       source_grantemp(iphs) = source_grantemp(iphs) - half*(zetan&
            & - abszetan)*fluctvmag2
       dissip_grantemp(iphs) = dissip_grantemp(iphs) + half*(zetan&
            & + abszetan)*fluctvmag2

       source_pres_grantemp(iphs) = source_pres_grantemp(iphs) - half*(zetanpres&
            & - abszetanpres)*fluctvmag2
       dissip_pres_grantemp(iphs) = dissip_grantemp(iphs) + half*(zetanpres&
            & + abszetanpres)*fluctvmag2

       source_visc_grantemp(iphs) = source_visc_grantemp(iphs) - half*(zetanvisc&
            & - abszetanvisc)*fluctvmag2
       dissip_visc_grantemp(iphs) = dissip_visc_grantemp(iphs) + half*(zetanvisc&
            & + abszetanvisc)*fluctvmag2
    end do

    source_grantemp(iphs) = two/three*source_grantemp(iphs)&
         &/real(npart, prcn)
    dissip_grantemp(iphs) = two/three*dissip_grantemp(iphs)&
         &/real(npart, prcn)

    source_pres_grantemp(iphs) = two/three*source_pres_grantemp(iphs)&
         &/real(npart, prcn)
    dissip_pres_grantemp(iphs) = two/three*dissip_pres_grantemp(iphs)&
         &/real(npart, prcn)
    
    source_visc_grantemp(iphs) = two/three*source_visc_grantemp(iphs)&
         &/real(npart, prcn)
    dissip_visc_grantemp(iphs) = two/three*dissip_visc_grantemp(iphs)&
         &/real(npart, prcn)

!!$    Write(*,*)'Source = ', source_grantemp(iphs), 'Source pres +  visc&
!!$         & = ',(source_pres_grantemp(iphs) +&
!!$         & source_visc_grantemp(iphs)-dissip_pres_grantemp(iphs) -&
!!$         & dissip_visc_grantemp(iphs))
!!$
!!$    Write(*,*)'Dissip = ', dissip_grantemp(iphs), 'Dissip pres +  visc&
!!$         & = ',dissip_pres_grantemp(iphs) + dissip_visc_grantemp(iphs)


  END SUBROUTINE compute_grantemp_source_dissipation

  SUBROUTINE compute_measvol_stats
        USE constants  
    IMPLICIT NONE
    CHARACTER(LEN=80) :: FILENAME, NUMVOLFNAME, TEMPFNAME, table_name
    INTEGER :: m, timestepcount, measvolcount, ICELLS, JCELLS, KCELLS&
         &, numvolfluctunit, i, j, k, idim, miscount, imeasvol,&
         & numvollineunit, tunit
    REAL(prcn) :: Lmeas, Lmeasmin, Lmeasmax, dmeasvol, charvol,&
         & charvolsq, Lmeaslin
    REAL(prcn) :: CONFIN_NUMBER, confin_numbersq,&
         & confin_volume, confin_volumesq, confin_velpar,&
         & confin_velper
    
    INTEGER :: nrbins, timestepmax, gofravgunit, ibin
    
    PARAMETER(nrbins = 200, timestepmax=100)
    
    REAL(prcn) :: gofr(nrbins), gofr_avg(nrbins), gofr_std(nrbins),&
         & rad_bin(nrbins), Ln, homog_number, homog_numberp, homog_numberm, num_analy
    
    LOGICAL :: filexist, isopen
    CHARACTER :: PER_CONF*8
    
!!$    Lmeasmin = 5.d0*dx
!!$    Lmeasmax = DOML(2)
!!$    dmeasvol = (Lmeasmax - Lmeasmin)/real(nmeasvols-1, prcn)
!!$    
!!$    do m = 1, nbody
!!$       do idim = 1, ndim
!!$          XC_PHYS(m,idim) = (XC(m,idim)-one)*dx
!!$       end do
!!$       RAD_PHYS(m) = RADBDY(m)*dx
!!$    end do
!!$
!!$    do imeasvol = 1, nmeasvols
!!$       Write(*,*)'measurement volume = ', imeasvol
!!$       Lmeaslin = Lmeasmin + (imeasvol-1) * dmeasvol
!!$       !Write(*,*)'Input measurement length = ', Lmeaslin
!!$       
!!$       !CALL compute_num_vol_in_measvol(Lmeaslin, imeasvol)
!!$       CALL compute_acc_phi_scatter(Lmeaslin,imeasvol)
!!$       !Lmeaslin = Lmeaslin/two
!!$    end do
!!$
!!$    STOP

    FILENAME = TRIM(RUN_NAME)//'_part_info'//'.rst'
    
    charvol = real(nbody, prcn)*pi*(dia_phys**3.d0)/6
    charvolsq = charvol**2.d0
    numdens = REAL(nbody,prcn)/(DOML(1)*DOML(2)*DOML(3))
    Num_analy = numdens * voldom

    Lmeasmin = dia_phys
    Lmeas = Lmeasmin
    Ln = (one/numdens)**(one/three)

    measvolcount = 0

    numvollineunit = getnewunit(minunitno, maxunitno)

    OPEN(numvollineunit, FILE=TRIM(RUN_NAME)//'_NUMVOLLINEPLOT.dat',&
         & form='formatted')

    gofravgunit = getnewunit(minunitno, maxunitno)
    
    OPEN(gofravgunit, FILE=TRIM(RUN_NAME)//'_GOFRAVG.dat',&
         & form='formatted')


    Lmeasmax = DOML(2)
    dmeasvol = (Lmeasmax - Lmeasmin)/real(nmeasvols-1, prcn)

    
    NUMBER_IN_MEASVOL = zero
    NUMBER_STDDEV = zero

    VEL_PAR_MEASVOL = zero
    VEL_PER_MEASVOL = zero

    VEL_PAR_STDDEV = zero
    VEL_PER_STDDEV = zero

    NUMBERSQ_IN_MEASVOL = zero

    VOLUME_IN_MEASVOL = zero
    VOLUME_STDDEV = zero

    VOLUMESQ_IN_MEASVOL = zero

    NUMBERSQ_STDDEV = zero
    VOLUMESQ_STDDEV = zero
    
    gofr_avg = zero
    gofr_std = zero

    numvolfluctunit = getnewunit(minunitno, maxunitno)
    
    timestepcount = 0
    miscount = 0
    
    OPEN(1050, FILE = FILENAME, STATUS='old', form='unformatted')
    Write(*,*)'READING THE FILE FOR GEOMETRIC STATISTICS ON MEAS VO&
         &L ', measvolcount , ' OF LENGTH = ', Lmeas
    
1250 continue    
    timestepcount = timestepcount+1
    READ(1050,END=2250)post_proc_time
    READ(1050,END=2250)xc(1:nbody, 1:3)
    READ(1050,END=2250)velbdy(1:nbody, 1:3)
    READ(1050,END=2250)force(1:nbody, 1:3)
    READ(1050,END=2250) pres(1:nbody, 1:3)
    READ(1050,END=2250) visc(1:nbody, 1:3)
!!$    READ(1050,END=2250) frame_vel(1:3)
!!$    READ(1050,END=2250) frame_accln(1:3)
!!$    READ(1050,END=2250) ufmean(1:3)
    !Write(*,*)'READING ', timestepcount
    if(MOD(timestepcount,500).eq.0)then
       miscount = miscount + 1
       CALL grid_nodes_insphere
       !Convert Positions and radii to physical units
       do m = 1, nbody
          do idim = 1, ndim
             XC_PHYS(m,idim) = (XC(m,idim)-one)*dx
          end do
          RAD_PHYS(m) = RADBDY(m)*dx
       end do
       
       Write(*,*)'COMPUTING NUMBER AND VOL FLUCTUATIONS FOR MIS =&
            & ', miscount
       
       do imeasvol = 1, nmeasvols
          Write(*,*)'measurement volume = ', imeasvol
          Lmeaslin = Lmeasmin + (imeasvol-1) * dmeasvol
          !Write(*,*)'Input measurement length = ', Lmeaslin
          
          !CALL compute_num_vol_in_measvol(Lmeaslin, imeasvol)
          CALL compute_acc_phi_scatter(Lmeaslin,imeasvol)
          !Lmeaslin = Lmeaslin/two
       end do
       CALL calculate_gofr_homog(nbody,xc(1:nbody,1:3), my, mx-1,&
            & nrbins, .TRUE., gofr(1:nrbins), rad_bin(1:nrbins))
    end if

    do ibin = 1, nrbins
       gofr_avg(ibin) = gofr_avg(ibin) + gofr(ibin)
       gofr_std(ibin) = gofr_std(ibin) + gofr(ibin)**2.d0
    end do
    
    goto 1250
    
2250 continue
    
    CLOSE(1050, status='keep')
    
    CALL GET_CONFIN(miscount,confin)

    do imeasvol = 1, nmeasvols
       NUMBER_IN_MEASVOL(imeasvol) =&
            & NUMBER_IN_MEASVOL(imeasvol)/real(miscount,prcn)
       
       NUMBER_STDDEV(imeasvol) = NUMBER_STDDEV(imeasvol)&
            &/real(miscount,prcn)
       
       confin_number = NUMBER_STDDEV(imeasvol) -&
            & (NUMBER_IN_MEASVOL(imeasvol))**2.d0
       
       confin_number = confin*DSQRT(confin_number/real(miscount&
            &,prcn))
       
       VEL_PAR_MEASVOL(imeasvol) =&
            & VEL_PAR_MEASVOL(imeasvol)/real(miscount,prcn)
       
       VEL_PER_MEASVOL(imeasvol) =&
            & VEL_PER_MEASVOL(imeasvol)/real(miscount,prcn)
       
       VEL_PAR_STDDEV(imeasvol) = VEL_PAR_STDDEV(imeasvol)&
            &/real(miscount,prcn)
       
       VEL_PER_STDDEV(imeasvol) = VEL_PER_STDDEV(imeasvol)&
            &/real(miscount,prcn)
       
       confin_velpar = VEL_PAR_STDDEV(imeasvol) -&
            & (VEL_PAR_MEASVOL(imeasvol))**2.d0
       
       confin_velpar = confin*DSQRT(confin_velpar/real(miscount&
            &,prcn))
       
       confin_velper = VEL_PER_STDDEV(imeasvol) -&
            & (VEL_PER_MEASVOL(imeasvol))**2.d0
       
       confin_velper = confin*DSQRT(confin_velper/real(miscount&
            &,prcn))
       
       
       NUMBERSQ_IN_MEASVOL(imeasvol) =&
            & NUMBERSQ_IN_MEASVOL(imeasvol)/real(miscount,prcn)
       
       NUMBERSQ_STDDEV(imeasvol) = NUMBERSQ_STDDEV(imeasvol)&
            &/real(miscount,prcn)
       
       NUMBERSQ_STDDEV(imeasvol) = NUMBERSQ_STDDEV(imeasvol) -&
            & NUMBERSQ_IN_MEASVOL(imeasvol)**2.d0
       
       confin_numbersq = confin*DSQRT(NUMBERSQ_STDDEV(imeasvol)&
            &/real(miscount,prcn))
       
       VOLUME_IN_MEASVOL(imeasvol) =&
            & VOLUME_IN_MEASVOL(imeasvol)/real(miscount,prcn)
       
       
       VOLUME_STDDEV(imeasvol) = VOLUME_STDDEV(imeasvol)&
            &/real(miscount,prcn)
       
       confin_volume = VOLUME_STDDEV(imeasvol) -&
            & (VOLUME_IN_MEASVOL(imeasvol))**2.d0
       
       confin_volume = confin*DSQRT(confin_volume/real(miscount&
            &,prcn))
       
       VOLUMESQ_IN_MEASVOL(imeasvol) =&
            & VOLUMESQ_IN_MEASVOL(imeasvol)/real(miscount,prcn)
       
       VOLUMESQ_STDDEV(imeasvol) = VOLUMESQ_STDDEV(imeasvol)&
                  &/real(miscount,prcn)
       VOLUMESQ_STDDEV(imeasvol) = VOLUMESQ_STDDEV(imeasvol) -&
            & VOLUMESQ_IN_MEASVOL(imeasvol)**2.d0
       confin_volumesq = confin*DSQRT(VOLUMESQ_STDDEV(imeasvol))&
            &/real(miscount,prcn)
       
       Lmeaslin = Lmeasmin + (imeasvol-1) * dmeasvol
       
       homog_number = numdens*(Lmeaslin**3.d0)
       homog_numberp = homog_number + DSQRT(homog_number)
       homog_numberm = homog_number - DSQRT(homog_number)
       Write(numvollineunit,'(20(2x,g17.5))')Lmeaslin**3.d0/voldom,&
            Lmeaslin, Lmeaslin/Ln,&
            & NUMBER_IN_MEASVOL(imeasvol)&
            &,confin_number, &
            & homog_number/num_analy, &
            & homog_numberp/num_analy, &
            & homog_numberm/num_analy, &
            & NUMBERSQ_IN_MEASVOL(imeasvol),&
            & confin_numbersq,&
            & VOLUME_IN_MEASVOL(imeasvol),&
            & confin_volume,&
            & VOLUMESQ_IN_MEASVOL(imeasvol),&
            & confin_volumesq,&
            & VEL_PAR_MEASVOL(imeasvol), confin_velpar,&
            & VEL_PER_MEASVOL(imeasvol), confin_velper
       
!!$             Write(*,*) ' Mean Number in meas vol = ',&
!!$                  & number_in_measvol(imeasvol), numdens*Lmeaslin&
!!$                  &**3.d0
    end do
    CLOSE(numvollineunit, status='keep')

    do ibin=1,nrbins
       gofr_avg(ibin) = gofr_avg(ibin)/real(miscount,prcn)
       gofr_std(ibin) = gofr_std(ibin)/real(miscount,prcn) -&
            & gofr_avg(ibin)**2.d0
       
       write(gofravgunit,'(10(E20.10,1x))') rad_bin(ibin)&
            &,rad_bin(ibin)*lybyd,gofr_avg(ibin),confin&
            &*DSQRT(gofr_std(ibin))/real(miscount,prcn)
       
    end do
    close(gofravgunit, status = "keep") 
 
  END SUBROUTINE compute_measvol_stats

  SUBROUTINE compute_auto_correlation
    IMPLICIT NONE
!!$             if(.not.compute_auto_corr)compute_auto_corr = grant_converged
!!$             compute_auto_corr = .TRUE.
!!$             if(compute_auto_corr)then
!!$                if(.not.base_state_saved)then
!!$                   if(ALLOCATED(velbdy0))then
!!$                      DEALLOCATE(velbdy0)
!!$                   end if
!!$                   ALLOCATE(velbdy0(nbody,ndim))
!!$                   do idim = 1, ndim
!!$                      mean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
!!$                   end do
!!$                   
!!$                   do idim = 1, ndim
!!$                      do m = 1, nbody 
!!$                         velbdy0(m,idim) = velbdy(m,idim) - mean_vel(idim)
!!$                      end do
!!$                   end do
!!$                   base_state_saved = .TRUE.
!!$                   
!!$                   do idim = 1, ndim
!!$                      rhoofs0(idim) = zero
!!$                      do m = 1, nbody
!!$                         rhoofs0(idim) = rhoofs0(idim) + velbdy0(m,idim)*velbdy0(m,idim) !(velbdy(m,idim)-mean_vel(idim))
!!$                      end do
!!$                      rhoofs0(idim) = rhoofs0(idim)/real(nbody,prcn)
!!$                   end do
!!$                   
!!$                   AUTO_CORR_SEP_TIME = ZERO
!!$                   WRITE(*,*)'MEAN VEL BASE = ', MEAN_VEL(1:ndim)
!!$                else
!!$                   AUTO_CORR_SEP_TIME = AUTO_CORR_SEP_TIME + dt
!!$                end if
!!$
!!$                do idim = 1, ndim
!!$                   rhoofs(idim) = zero
!!$                   lag_str_func(idim) = zero
!!$                   mean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
!!$                   do m = 1, nbody
!!$                      rhoofs(idim) = rhoofs(idim) + velbdy0(m,idim)*(velbdy(m,idim)-mean_vel(idim))
!!$                      lag_str_func(idim) = lag_str_func(idim) + ((velbdy(m,idim)-mean_vel(idim))-velbdy0(m,idim))**2.d0
!!$                   end do
!!$                   rhoofs(idim) = rhoofs(idim)/real(nbody,prcn)
!!$                   lag_str_func(idim) = lag_str_func(idim)/real(nbody,prcn)
!!$                   !rhoofs(idim) = rhoofs(idim)/rhoofs0(idim)
!!$                end do
!!$                WRITE(unit_auto_corr, '(10(2x,g17.8))')
    !! auto_corr_sep_time, auto_corr_sep_time/t_conv,
    !! SUM(rhoofs(1:ndim))/SUM(rhoofs0(1:ndim)),
    !! lag_str_func(1:ndim), (SUM(lag_str_func(1:ndim)))
    !!/(SUM(rhoofs0(1:ndim))) 
!!$          end if


  END SUBROUTINE compute_auto_correlation


  SUBROUTINE compute_moment_of_inertia_principal_coordinates
#if 0
    IMPLICIT NONE
    INTEGER :: m, idim, iphs, jdim, njacobirot
    REAL(prcn) :: center_of_mass(ndim), totmass, I_Cart(ndim,ndim)

    center_of_mass = zero
    totmass = zero
    do m = 1, nbody
       iphs = part_array(m)%iphs
       totmass = totmass + phase_mass(iphs)
       !totmass = totmass + 1.0
       do idim = 1, ndim
          center_of_mass(idim) = center_of_mass(idim) + phase_mass(iphs)&
               &*xc(m,idim)
!!$          center_of_mass(idim) = center_of_mass(idim) + &
!!$               &xc(m,idim)
       end do
    end do
    do idim = 1, ndim
       center_of_mass(idim) = center_of_mass(idim)/totmass
    end do
    !    Write(*,*) ' Center of mass = ', center_of_mass(1:ndim)
    do idim = 1, ndim
       do jdim = idim, ndim
          I_Cart(idim,jdim) = zero
          do m = 1, nbody
             iphs = part_array(m)%iphs
             I_Cart(idim,jdim) = I_Cart(idim,jdim) +&
                  & phase_mass(iphs)*(xc(m,idim)&
                  &-center_of_mass(idim))*(xc(m,jdim)&
                  &-center_of_mass(jdim))
!!$             I_Cart(idim,jdim) = I_Cart(idim,jdim) +&
!!$                  & (xc(m,idim)&
!!$                  &-center_of_mass(idim))*(xc(m,jdim)&
!!$                  &-center_of_mass(jdim))
          end do
       end do
    end do
    !Moment of inertia is symmetric
    do jdim = 1, ndim
       do idim = 1, jdim-1
          I_Cart(jdim,idim) = I_Cart(idim,jdim)
       end do
    end do

    !Now compute the principal coordinates and eigen values of this
    ! moment of inertia matrix
    CALL jacobi(I_Cart(1:ndim,1:ndim),ndim,ndim,eigenvalue(1:ndim)&
         &,eigen_vectors(1:ndim,1:ndim),njacobirot)
#endif
  END SUBROUTINE compute_moment_of_inertia_principal_coordinates



  SUBROUTINE compute_geometric_statistics
    USE constants  
    IMPLICIT NONE
    CHARACTER(LEN=80) :: FILENAME, NUMVOLFNAME, TEMPFNAME, table_name
    INTEGER :: m, timestepcount, measvolcount, ICELLS, JCELLS, KCELLS&
         &, numvolfluctunit, i, j, k, idim, miscount, imeasvol,&
         & numvollineunit, tunit
    REAL(prcn) :: Lmeas, Lmeasmin, Lmeasmax, dmeasvol, charvol,&
         & charvolsq, Lmeaslin
    REAL(prcn) :: CONFIN_NUMBER, confin_numbersq,&
         & confin_volume, confin_volumesq, confin_velpar,&
         & confin_velper
    
    INTEGER :: nrbins, timestepmax, gofravgunit, ibin
    
    PARAMETER(nrbins = 200, timestepmax=100)
    
    REAL(prcn) :: gofr(nrbins), gofr_avg(nrbins), gofr_std(nrbins),&
         & rad_bin(nrbins), Ln, homog_number, homog_numberp, homog_numberm, num_analy
    
    LOGICAL :: ONE_MEAS_VOL_DONE,filexist, isopen
    CHARACTER :: PER_CONF*8
    
    FILENAME = TRIM(RUN_NAME)//'_part_info'//'.rst'
    
    charvol = real(nbody, prcn)*pi*(dia_phys**3.d0)/6
    charvolsq = charvol**2.d0
    numdens = REAL(nbody,prcn)/(DOML(1)*DOML(2)*DOML(3))
    Num_analy = numdens * voldom

    Lmeasmin = 5.d0*dx
    Lmeas = Lmeasmin
    Ln = (one/numdens)**(one/three)

    measvolcount = 0

    

    Lmeasmax = DOML(2)
    dmeasvol = (Lmeasmax - Lmeasmin)/real(nmeasvols-1, prcn)

    numvolfluctunit = getnewunit(minunitno, maxunitno)
    Do While(Lmeas.ge.Lmeasmin)
       
       measvolcount = measvolcount + 1
       
       Write(TEMPFNAME, FMT="(I2.2)")measvolcount
       NUMVOLFNAME = TRIM(RUN_NAME)//'_NUMVOLFLUCTS_MEASVOL'&
            &//TRIM(TEMPFNAME)//'.dat'

       timestepcount = 0
       miscount = 0
       OPEN(unit = numvolfluctunit, FILE=TRIM(NUMVOLFNAME), FORM='forma&
            &tted')

       IF(ASSOCIATED(COARSE_GRID%XE))DEALLOCATE(COARSE_GRID%XE)
       IF(ASSOCIATED(COARSE_GRID%YN))DEALLOCATE(COARSE_GRID%YN)
       IF(ASSOCIATED(COARSE_GRID%ZT))DEALLOCATE(COARSE_GRID%ZT)

       CALL INITIALIZE_GRID(Lmeas, ICELLS, JCELLS, KCELLS)

       IF(ASSOCIATED(COARSE_GRID%meannum))DEALLOCATE(COARSE_GRID&
            &%meannum)
       IF(ASSOCIATED(COARSE_GRID%meannumsq))DEALLOCATE(COARSE_GRID&
            &%meannumsq)

       IF(ASSOCIATED(COARSE_GRID%meanvol))DEALLOCATE(COARSE_GRID&
            &%meanvol)
       IF(ASSOCIATED(COARSE_GRID%meanvolsq))DEALLOCATE(COARSE_GRID&
            &%meanvolsq)

       ALLOCATE(COARSE_GRID%meannum(ICELLS, JCELLS, KCELLS),&
            & COARSE_GRID%meanvol(ICELLS, JCELLS, KCELLS))
       ALLOCATE(COARSE_GRID%meannumsq(ICELLS, JCELLS, KCELLS),&
            & COARSE_GRID%meanvolsq(ICELLS, JCELLS, KCELLS))

       Do k = 1, kcells
          Do j = 1, jcells
             Do i = 1, icells
                COARSE_GRID%meannum(i,j,k) = zero
                COARSE_GRID%meanvol(i,j,k) = zero
                COARSE_GRID%meannumsq(i,j,k) = zero
                COARSE_GRID%meanvolsq(i,j,k) = zero
             End Do
          End Do
       End Do

       OPEN(1050, FILE = FILENAME, STATUS='old', form='unformatted')
       Write(*,*)'READING THE FILE FOR GEOMETRIC STATISTICS ON MEAS VO&
            &L ', measvolcount , ' OF LENGTH = ', Lmeas

1250   continue    
       timestepcount = timestepcount+1
       READ(1050,END=2250)post_proc_time
       READ(1050,END=2250)xc(1:nbody, 1:3)
       READ(1050,END=2250)velbdy(1:nbody, 1:3)
       READ(1050,END=2250)force(1:nbody, 1:3)
       READ(1050,END=2250) pres(1:nbody, 1:3)
       READ(1050,END=2250) visc(1:nbody, 1:3)
       READ(1050,END=2250) frame_vel(1:3)
       READ(1050,END=2250) frame_accln(1:3)
       READ(1050,END=2250) ufmean(1:3)
       !Write(*,*)'READING ', timestepcount
       if(MOD(timestepcount,100).eq.0)then
          miscount = miscount + 1
          CALL grid_nodes_insphere
          !Convert Positions and radii to physical units
          do m = 1, nbody
             do idim = 1, ndim
                XC_PHYS(m,idim) = (XC(m,idim)-one)*dx
             end do
             RAD_PHYS(m) = RADBDY(m)*dx
          end do

          Write(*,*)'COMPUTING NUMBER AND VOL FLUCTUATIONS FOR MIS =&
               & ', miscount
          !CALL compute_number_volume_fluct_fields(Lmeas)
          !STOP
       end if
       goto 1250

2250   continue
       
       CLOSE(1050, status='keep')
       write(numvolfluctunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "NUM" ', ' "NUMSQ" &
            &',' "VOL" ',' "VOLSQ" ' 
       write(numvolfluctunit,*)'ZONE F=POINT, I=', ICELLS,  ', J=', JCELLS, ', K&
            &=', KCELLS
       Do k = 1, kcells
          Do j = 1, jcells
             Do i = 1, icells
                COARSE_GRID%meannum(i,j,k) = COARSE_GRID%meannum(i,j&
                     &,k)/real(miscount, prcn)
                COARSE_GRID%meanvol(i,j,k) = COARSE_GRID%meanvol(i,j&
                     &,k)/real(miscount, prcn)

                COARSE_GRID%meannumsq(i,j,k) = COARSE_GRID%meannumsq(i,j&
                     &,k)/real(miscount, prcn)
                COARSE_GRID%meanvolsq(i,j,k) = COARSE_GRID%meanvolsq(i,j&
                     &,k)/real(miscount, prcn)

                Write(numvolfluctunit, *)real(i), real(j), real(k),&
                     & COARSE_GRID%meannum(i,j,k), COARSE_GRID&
                     &%meannumsq(i,j,k), COARSE_GRID%meanvol(i,j&
                     &,k)/charvol, (COARSE_GRID%meanvolsq(i,j&
                     &,k))/charvolsq
             End Do
          End Do
       End Do

       CLOSE(numvolfluctunit, status='keep')

       Write(*,*)'Total Mean number = ', SUM(COARSE_GRID%meannum(:,:,:))
       Write(*,*)'Total Mean volume = ', SUM(COARSE_GRID%meanvol(:,:,:))
       
       Lmeas = Lmeas/two
    End Do

  END SUBROUTINE compute_geometric_statistics

  SUBROUTINE initialize_grid(Lmeas, ICELLS, JCELLS, KCELLS)
    IMPLICIT NONE
    REAL(prcn), Intent(in) :: Lmeas
    INTEGER, INTENT(out) :: ICELLS, JCELLS, KCELLS
    INTEGER :: cx, cy, cz, i, imin, jmin,&
         & kmin, ifinew,ifinee, jfines, jfinen, kfineb,&
         & kfinet, j, k

    cx = MAX(INT(DOML(1)/Lmeas)+1,2)
    cy = MAX(INT(DOML(2)/Lmeas)+1,2)
    cz = MAX(INT(DOML(3)/Lmeas)+1,2)

    COARSE_GRID%dx = DOML(1)/REAL(cx-1, prcn)
    COARSE_GRID%dy = DOML(2)/REAL(cy-1, prcn)
    COARSE_GRID%dz = DOML(3)/REAL(cz-1, prcn)

    COARSE_GRID%cx = cx
    COARSE_GRID%cy = cy
    COARSE_GRID%cz = cz

    ALLOCATE(coarse_grid%XE(cx-1),coarse_grid%YN(cy-1),coarse_grid&
         &%ZT(cz-1))

    ! ICELLS,JCELLS,KCELLS are the indices of the last physical CELL
    ICELLS = CX-1
    JCELLS = CY-1
    KCELLS = CZ-1

    ! IMIN,JMIN,KMIN are the indices of the first physical CELL
    IMIN = 1
    JMIN = 1
    KMIN = 1

    COARSE_GRID%XE(IMIN) = COARSE_GRID%dx
    COARSE_GRID%YN(JMIN) = COARSE_GRID%dy
    COARSE_GRID%ZT(KMIN) = COARSE_GRID%dz
    do i = IMIN+1, ICELLS
       COARSE_GRID%XE(I) = COARSE_GRID%XE(I-1) + COARSE_GRID%dx
    end do
    do i = JMIN+1, JCELLS
       COARSE_GRID%YN(I) = COARSE_GRID%YN(I-1) + COARSE_GRID%dy
    end do
    do i = KMIN+1, KCELLS
       COARSE_GRID%ZT(I) = COARSE_GRID%ZT(I-1) + COARSE_GRID%dz
    end do
    
  END SUBROUTINE initialize_grid
  


  SUBROUTINE compute_num_vol_in_measvol(Lmeaslin, imeasvol)
    USE constants  
    IMPLICIT NONE
    INTEGER, Intent(in) :: imeasvol
    REAL(prcn), Intent(in) :: Lmeaslin
    REAL(prcn) :: XW, XE, YS, YN, ZB, ZT, temp_pos(ndim), volume
    INTEGER :: m, number, idim
    INTEGER :: ifinew,ifinee, jfines, jfinen, kfineb,&
         & kfinet, ii, jj, kk, PIJK(ndim)
    LOGICAL :: XINMEASVOL, YINMEASVOL, ZINMEASVOL
    REAL(prcn) :: num_analy, alpha2, ssfm, homog_num, ssfm_vol, homog_vol, vol_anal, vel_measvol(ndim)

    Num_analy = numdens * voldom
    homog_num = numdens * (Lmeaslin**3.d0)

    !homog_vol = homog_num*pi*dia_phys**3.d0/6.d0
    homog_vol = homog_num*pi*(dbydx**3.d0)/6.d0

    XW = DOML(1)/two - Lmeaslin/two
    XE = DOML(1)/two + Lmeaslin/two

    YS = DOML(2)/two - Lmeaslin/two
    YN = DOML(2)/two + Lmeaslin/two

    ZB = DOML(3)/two - Lmeaslin/two
    ZT = DOML(3)/two + Lmeaslin/two

    number = 0
    vel_measvol = zero
    
    do M = 1, nbody
       do idim = 1, ndim
          temp_pos(idim) = XC_PHYS(m,idim)
       end do
       XINMEASVOL = (temp_pos(1).ge.XW).and.(temp_pos(1).le.XE)
       YINMEASVOL = (temp_pos(2).ge.YS).and.(temp_pos(2).le.YN)
       ZINMEASVOL = (temp_pos(3).ge.ZB).and.(temp_pos(3).le.ZT)
       if(XINMEASVOL.AND.YINMEASVOL.AND.ZINMEASVOL)then
          number = number + 1
          vel_measvol(1:ndim) = vel_measvol(1:ndim) + velbdy(m,1:ndim)
       end if
    end do


    KFINEB = INT(ZB/dz) + 1
    KFINET = INT(ZT/dz) + 1
    JFINES = INT(YS/dy) + 1
    JFINEN = INT(YN/dy) + 1
    IFINEW = INT(XW/dx) + 1
    IFINEE = INT(XE/dx) + 1

!!$
!!$    KFINEB = CEILING(ZB/dz)!INT(ZB/dz) + 1
!!$    KFINET = FLOOR(ZT/dz) !INT(ZT/dz) + 1
!!$    JFINES = CEILING(YS/dy) !INT(YS/dy) + 1
!!$    JFINEN = FLOOR(YN/dy)!INT(YN/dy) + 1
!!$    IFINEW = CEILING(XW/dx) !INT(XW/dx) + 1
!!$    IFINEE = FLOOR(XE/dx) !INT(XE/dx) + 1
    
    volume = zero
    
    do KK = KFINEB, KFINET
       do JJ = JFINES, JFINEN
          do II = IFINEW, IFINEE
             if((KK.le.mz).and.(JJ.le.my).and.(II.le.mx1))then
                if(.not.fluid_atijk(II,JJ,KK))then
                   !volume = volume + dx**3.d0
                   volume = volume + one**3.d0
                end if
             end if
          end do
       end do
    end do
    
    
    NUMBER_IN_MEASVOL(imeasvol) = NUMBER_IN_MEASVOL(imeasvol)+&
         & real(number,prcn)/num_analy
    
    NUMBER_STDDEV(imeasvol) = NUMBER_STDDEV(imeasvol) + (real(number&
         &,prcn)/num_analy)**2.d0

    if(number.gt.0)then
       vel_measvol(1:ndim) = vel_measvol(1:ndim)/real(number,prcn)
    end if

    VEL_PAR_MEASVOL(imeasvol) = vel_measvol(1) 
    VEL_PER_MEASVOL(imeasvol) = vel_measvol(2)

    VEL_PAR_STDDEV(imeasvol) = VEL_PAR_STDDEV(imeasvol) + vel_measvol(1)**2.d0
    VEL_PER_STDDEV(imeasvol) = VEL_PAR_STDDEV(imeasvol) + vel_measvol(2)**2.d0

    ssfm = (number**2.d0 - number)/(homog_num)**2.d0
    
    ssfm_vol = (volume**2.d0 - volume)/(homog_vol)**2.d0

#if 1
    if(imeasvol.eq.1) then
       Write(*,*) ' Number in meas vol = ', number, num_analy
       Write(*,*) ' SSFM = ', ssfm
       Write(*,*) ' Volume in meas vol = ', volume, homog_vol
       Write(*,*) ' SSFM = ', ssfm_vol

    end if
#endif

    
    NUMBERSQ_IN_MEASVOL(imeasvol) = NUMBERSQ_IN_MEASVOL(imeasvol)+ ssfm
    
    NUMBERSQ_STDDEV(imeasvol) = NUMBERSQ_STDDEV(imeasvol) + ssfm**2.d0
    
    VOLUME_IN_MEASVOL(imeasvol) = VOLUME_IN_MEASVOL(imeasvol)+ volume&
         &/homog_vol

    VOLUME_STDDEV(imeasvol) = VOLUME_STDDEV(imeasvol)+ (volume&
         &/homog_vol)**2.d0
    
    VOLUMESQ_IN_MEASVOL(imeasvol) = VOLUMESQ_IN_MEASVOL(imeasvol)+&
         & ssfm_vol
    VOLUMESQ_STDDEV(imeasvol) = VOLUMESQ_STDDEV(imeasvol) + ssfm_vol&
         &**2.d0
    
    !if(imeasvol.eq.nmeasvols) Write(*,*) ' TOtal Number in meas vol = ',&
    !             & number_in_measvol(imeasvol)

  END SUBROUTINE compute_num_vol_in_measvol

  SUBROUTINE compute_acc_phi_scatter(Lmeaslin, imeasvol)
    USE constants 
    USE dependent_functions
    USE nlmainarrays, velr=>ubcp
    IMPLICIT NONE
    INTEGER, Intent(in) :: imeasvol
    REAL(prcn), Intent(in) :: Lmeaslin
    REAL(prcn) :: XW, XE, YS, YN, ZB, ZT, temp_pos(ndim), volume
    REAL(prcn) :: part_vol(nbody)
    INTEGER :: m, number, idim
    INTEGER :: ifinew,ifinee, jfines, jfinen, kfineb,&
         & kfinet, ii, jj, kk, PIJK(ndim)
    LOGICAL :: XINMEASVOL, YINMEASVOL, ZINMEASVOL
    REAL(prcn) :: num_analy, alpha2, ssfm, homog_num, ssfm_vol,&
         & homog_vol, vol_anal, vel_measvol(ndim), xlr(ndim) ,&
         & xll(ndim), acc_measvol(ndim), fluid_vel(ndim),&
         & volume_fraction, measvolslip(ndim), measvolslipmod,&
         & Rem_measvol, volume_part, F_IBM, measvol_drag,&
         & force_drag_law
    INTEGER :: cor_min(ndim), cor_max(ndim), imin, imax, jmin, jmax,&
         & kmin, kmax, i, j, k, measvolcountfl
    
    Num_analy = numdens * voldom
    homog_num = numdens * (Lmeaslin**3.d0)
    
    !homog_vol = homog_num*pi*dia_phys**3.d0/6.d0
    homog_vol = homog_num*pi*(dbydx**3.d0)/6.d0
    
    XW = DOML(1)/two - Lmeaslin/two
    XE = DOML(1)/two + Lmeaslin/two
    
    YS = DOML(2)/two - Lmeaslin/two
    YN = DOML(2)/two + Lmeaslin/two
    
    ZB = DOML(3)/two - Lmeaslin/two
    ZT = DOML(3)/two + Lmeaslin/two
    
    number = 0
    vel_measvol = zero
    acc_measvol = zero
    
    Do m = 1, nbody
       part_vol(m) = zero
       do idim = 1, ndim
          xlr(idim) = xc(m,idim)  + radbdy(m)
          xll(idim) = xc(m,idim)  - radbdy(m) 
       end do
       do idim = 1, ndim 
          cor_min(idim)  = ceiling(xll(idim))
          cor_max(idim) = floor(xlr(idim)) 
       end do
       imin = cor_min(1)
       imax = cor_max(1)
       jmin = cor_min(2)
       jmax = cor_max(2)
       kmin = cor_min(3)
       kmax = cor_max(3)

       do i = imin, imax 
          ii = i
          if(i.lt.1.and.intx_per) ii = mxf+i-1
          if(i.gt.mxf-1.and.intx_per) ii = i-(mxf-1)
          do j = jmin, jmax 
             jj = j 
             if(j.lt.1.and.inty_per) jj = my+j
             if(j.gt.my.and.inty_per) jj = j-my
             do k = kmin, kmax 
                kk = k 
                if(k.lt.1.and.intz_per) kk = mz+k
                if(k.gt.mz.and.intz_per) kk = k-mz 
                
                if(fluid_atijk(ii,jj,kk))then
                   temp_pos(1) = (ii-1)*dx
                   temp_pos(2) = (jj-1)*dy
                   temp_pos(3) = (kk-1)*dz
                   XINMEASVOL = (temp_pos(1).ge.XW).and.(temp_pos(1).le.XE)
                   YINMEASVOL = (temp_pos(2).ge.YS).and.(temp_pos(2).le.YN)
                   ZINMEASVOL = (temp_pos(3).ge.ZB).and.(temp_pos(3).le.ZT)
                   if(XINMEASVOL.AND.YINMEASVOL.AND.ZINMEASVOL)then
                      part_vol(m) = part_vol(m) + dx**3.d0
                   end if
                end if
             end do
          end do
       end do
       
    end Do

    ! Calculating total volume of intersection in another way
    do M = 1, nbody
       do idim = 1, ndim
          temp_pos(idim) = XC_PHYS(m,idim)
       end do
       XINMEASVOL = (temp_pos(1).ge.XW).and.(temp_pos(1).le.XE)
       YINMEASVOL = (temp_pos(2).ge.YS).and.(temp_pos(2).le.YN)
       ZINMEASVOL = (temp_pos(3).ge.ZB).and.(temp_pos(3).le.ZT)
       if(XINMEASVOL.AND.YINMEASVOL.AND.ZINMEASVOL)then
          number = number + 1
          vel_measvol(1:ndim) = vel_measvol(1:ndim) + velbdy(m,1:ndim)
       end if
       volume_part = pi*(two*radbdy(m)*dx)**3.d0/6.d0
       acc_measvol(1:ndim) = acc_measvol(1:ndim) + part_vol(m)&
            &*(pres(m,1:ndim) + visc(m,1:ndim))/Lmeaslin**3.d0 !volume_part
    end do
    

    KFINEB = INT(ZB/dz) + 1
    KFINET = INT(ZT/dz) + 1
    JFINES = INT(YS/dy) + 1
    JFINEN = INT(YN/dy) + 1
    IFINEW = INT(XW/dx) + 1
    IFINEE = INT(XE/dx) + 1
    
    
    volume = zero
    measvolcountfl = 0
    fluid_vel(1:ndim) = zero
    
    CALL calc_velreal(velr) 
    
    do KK = KFINEB, KFINET
       do JJ = JFINES, JFINEN
          do II = IFINEW, IFINEE
             if((KK.le.mz).and.(JJ.le.my).and.(II.le.mx1))then
                if(.not.fluid_atijk(II,JJ,KK))then
                   volume = volume + dx**3.d0
                   !volume = volume + one**3.d0
                else
                   measvolcountfl = measvolcountfl + 1
                   fluid_vel(1:ndim) = fluid_vel(1:ndim) + velr(ii,jj,kk,1:ndim)
                end if
             end if
          end do
       end do
    end do
    
    
    volume_fraction = volume/(Lmeaslin**3.d0)
    
    if(measvolcountfl.gt.0)then
       fluid_vel(1:ndim) = fluid_vel(1:ndim)/real(measvolcountfl,prcn)
    end if
    
    if(number.gt.0)then
       vel_measvol(1:ndim) = vel_measvol(1:ndim)/real(number,prcn)
    end if
    
    measvolslip(1:ndim) = fluid_vel(1:ndim) - vel_measvol(1:ndim)
    measvolslipmod = DSQRT(DOT_PRODUCT(measvolslip(1:ndim), measvolslip(1:ndim)))
    
    Rem_measvol = measvolslipmod*(one-volume_fraction)/vis
    
    !Write(*,*)'MEASVOL : ', Lmeaslin, Rem_measvol

    CALL compute_ibm_drag(volume_fraction, Re, F_IBM)
    
    if(volume.gt.zero)then
       acc_measvol(1:ndim) = acc_measvol(1:ndim)/volume_fraction
    end if

    measvol_drag = DSQRT(DOT_PRODUCT(acc_measvol(1:ndim), acc_measvol(1:ndim)))
    measvol_drag = measvol_drag/(3.d0*pi*vis*ucharmod*dia_phys)
    
    force_drag_law = F_IBM
    
    VEL_PAR_MEASVOL(imeasvol) = vel_measvol(1) 
    VEL_PER_MEASVOL(imeasvol) = vel_measvol(2)

    VEL_PAR_STDDEV(imeasvol) = VEL_PAR_STDDEV(imeasvol) + vel_measvol(1)**2.d0
    VEL_PER_STDDEV(imeasvol) = VEL_PAR_STDDEV(imeasvol) + vel_measvol(2)**2.d0


#if 1
    if(imeasvol.eq.1) then
       Write(*,*) ' Number in meas vol = ', number, num_analy
       Write(*,*) ' SSFM = ', ssfm
       Write(*,*) ' Volume in meas vol = ', volume, homog_vol
       Write(*,*) ' SSFM = ', ssfm_vol

    end if
#endif

    Write(accphiunit, '(4(2x,g17.8))')volume_fraction, Rem_measvol, measvol_drag, force_drag_law     
    

  END SUBROUTINE compute_acc_phi_scatter

  SUBROUTINE compute_ibm_drag(phi, rem, F)
    IMPLICIT NONE
    Real(prcn), Intent(in) :: phi, rem
    Real(prcn), Intent(out) :: F
    Real(prcn) :: RE, FISOL, c, F0, F1

       RE =  Rem
       c = phi
       FISOL = 1.d0 + 0.15*(RE**0.687)
       FISOL = FISOL/(1-c)**3.d0
       F0 = 5.813*c/(1-c)**3.d0 + 0.485*c**(1.0/3.0)/(1-c)**4.d0
       F1 = RE*(c**3.d0)*(0.954 + 0.607*(c**3.d0/(1-c)**2.d0))
       F = FISOL + F0 + F1
   END SUBROUTINE compute_ibm_drag

  SUBROUTINE compute_number_volume_fluct_fields(Lmeas)
    USE constants  
    IMPLICIT NONE
    REAL(prcn), Intent(in) :: Lmeas
    INTEGER :: cx, cy, cz, i, imin, imax, jmin,&
         & jmax, kmin, kmax, ifinew,ifinee, jfines, jfinen, kfineb,&
         & kfinet, j, k, ii, jj, kk, PIJK(ndim), m, idim
    REAL(prcn) :: xe, xw, ys, yn, zb, zt, temp_pos(ndim)
    REAL(prcn), ALLOCATABLE, DIMENSION(:,:,:) :: volume
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: number
    LOGICAL :: ALREADY_COUNTED(mx1,my,mz)

    do k = 1, mz
       do j = 1, my
          do i = 1, mx1
             ALREADY_COUNTED(i,j,k) = .FALSE.
          end do
       end do
    end do

    cx = COARSE_GRID%cx 
    cy = COARSE_GRID%cy 
    cz = COARSE_GRID%cz 

    ! IMIN1,JMIN1,KMIN1 are the indices of the first physical CELL
    IMIN = 1
    JMIN = 1
    KMIN = 1

    ! IMAX1,JMAX1,KMAX1 are the indices of the last physical CELL
    IMAX = CX-1
    JMAX = CY-1
    KMAX = CZ-1

!!$    Write(*,*) 'IMAX = ', IMAX, JMAX, KMAX
!!$    Write(*,*) "COARSE GRID dx = ", COARSE_GRID%dx, COARSE_GRID%dx, COARSE_GRID%dx

    ALLOCATE(Number(IMAX,JMAX,KMAX), Volume(IMAX, JMAX, KMAX))

    Do k = KMIN, KMAX
       Do j = JMIN, JMAX
          Do i = IMIN, IMAX
             Number(i,j,k) = 0
             Volume(i,j,k) = zero
          End Do
       End Do
    End Do

    do M = 1, nbody
       do idim = 1, ndim
          temp_pos(idim) = XC_PHYS(m,idim)
       end do
       PIJK(1) = MIN(INT(temp_pos(1)/COARSE_GRID%dx)+1, IMAX)
       PIJK(2) = MIN(INT(temp_pos(2)/COARSE_GRID%dy)+1, JMAX)
       PIJK(3) = MIN(INT(temp_pos(3)/COARSE_GRID%dz)+1, KMAX)

       Number(PIJK(1), PIJK(2), PIJK(3)) = Number(PIJK(1), PIJK(2),&
            & PIJK(3)) + 1
    end do



!!$    Write(*,*)'NUMber =', Number(:,:,:)
    do k = KMIN, KMAX
       ZT = COARSE_GRID%ZT(K)
       ZB = ZT - COARSE_GRID%dz
       !Write(*,*)'ZT = ', ZT, ZB
       KFINEB = INT(ZB/dz) + 1
       KFINET = INT(ZT/dz) + 1
       !Write(*,*) ' KFINES = ', KFINEB, KFINET
       do j = JMIN, JMAX
          YN = COARSE_GRID%YN(J)
          YS = YN - COARSE_GRID%dy
          JFINES = INT(YS/dy) + 1
          JFINEN = INT(YN/dy) + 1
          !Write(*,*) ' JFINES = ', JFINES, JFINEN
          do i = IMIN, IMAX
             XE = COARSE_GRID%XE(I)
             XW = XE - COARSE_GRID%dx
             IFINEW = INT(XW/dx) + 1
             IFINEE = INT(XE/dx) + 1
             !Write(*,*) ' IFINES =', IFINEW, IFINEE
!!$             Write(*,*) ' JFINES =', JFINES, JFINEN
!!$             Write(*,*) ' KFINES =', KFINEB, KFINET
             do KK = KFINEB, KFINET
                do JJ = JFINES, JFINEN
                   do II = IFINEW, IFINEE
                      if((KK.le.mz).and.(JJ.le.my).and.(II.le.mx1))then
                         if(.not.fluid_atijk(II,JJ,KK))then
                            if(.not.ALREADY_COUNTED(II,JJ,KK))then
                               Volume(i,j,k) = Volume(i,j,k) + dx&
                                    &**3.d0
                               ALREADY_COUNTED(II,JJ,KK) = .TRUE.
                            end if
                         end if
                      end if
                   end do
                end do
             end do
             !Write(*,*) ' I = ', i, ' J = ', j, ' K = ', k, 'Vol = ',&
             !     & volume(i,j,k), ' Num = ', number(i,j,k)
             COARSE_GRID%meannum(i,j,k) = COARSE_GRID%meannum(i,j,k) &
                  &+ Number(i,j,k)
             COARSE_GRID%meannumsq(i,j,k) = COARSE_GRID%meannumsq(i,j,k) &
                  &+ Number(i,j,k)**2

             COARSE_GRID%meanvol(i,j,k) = COARSE_GRID%meanvol(i,j,k) &
                  &+ Volume(i,j,k)
             COARSE_GRID%meanvolsq(i,j,k) = COARSE_GRID%meanvolsq(i,j,k) &
                  &+ Volume(i,j,k)**2.d0

          end do
       end do
    end do
!!$    Write(*,*)'Volume =', volume(:,:,:)
!!$    Write(*,*)'Total number = ', SUM(number(:,:,:))
!!$    Write(*,*)'Total volume = ', SUM(volume(:,:,:))
!!$    Write(*,*)'Mean NUMber =', COARSE_GRID%meannum(:,:,:)
!!$    Write(*,*)'Mean Volume =', COARSE_GRID%meanvol(:,:,:)

    DEALLOCATE(Number, Volume)

  END SUBROUTINE compute_number_volume_fluct_fields

  SUBROUTINE compute_steady_state_statistics
    IMPLICIT NONE

    INTEGER :: iphs, phasespaceunit

    phasespaceunit = getnewunit(minunitno, maxunitno)
    OPEN(phasespaceunit, FILE =TRIM(RUN_NAME)//'_nondim_phase_space.da&
         &t' , STATUS='replace', form='formatted')


    CALL compute_mean_and_fluct_quantites

    CALL compute_granular_temperature(.FALSE.)

    do iphs = 1, nphases
       nondim_grantemp(iphs) = granular_temperature(iphs) !-&
       !& steady_granular_temp(iphs))/steady_granular_temp(iphs)

!!$       nondim_grantemp(iphs) = (granular_temperature(iphs) -&
!!$            & steady_granular_temp(iphs))/steady_granular_temp(iphs)

       CALL compute_grantemp_source_dissipation(iphs)       

       nondim_source_grantemp(iphs) = (source_grantemp(iphs)&
            &-steady_source_grantemp(iphs))&
            &/steady_source_grantemp(iphs)

       nondim_dissip_grantemp(iphs) = (dissip_grantemp(iphs)&
            &-steady_dissip_grantemp(iphs))/steady_dissip_grantemp(iphs)
    end do

    WRITE(phasespaceunit, '(9(2x&
         &,g17.5))')(nondim_grantemp(iphs),&
         & nondim_source_grantemp(iphs), nondim_dissip_grantemp(iphs), iphs = 1&
         &,nphases)




    CLOSE(phasespaceunit, status='keep')

  END SUBROUTINE compute_steady_state_statistics
  
  SUBROUTINE ALLOCATE_MOVING_PART_RELATED
    IMPLICIT NONE
    ALLOCATE(granular_temperature(nphases),char_force(nphases),&
         & meanforcemod(nphases), phase_mass(nphases), acclnstddev_mod(nphases),vivjzioft(nphases)&
         &,vivjetaoft(nphases),aivjzioft(nphases),&
         & aivjetaoft(nphases), source_grantemp(nphases),&
         & dissip_grantemp(nphases),&
         & steady_granular_temp(nphases),&
         & steady_source_grantemp(nphases),&
         & steady_dissip_grantemp(nphases),nondim_grantemp(nphases),&
         & nondim_source_grantemp(nphases)&
         &,nondim_dissip_grantemp(nphases))
    ALLOCATE(source_pres_grantemp(nphases),&
         & source_visc_grantemp(nphases),&
         & dissip_pres_grantemp(nphases),&
         & dissip_visc_grantemp(nphases), RAD_PHYS(nbody))
    ALLOCATE(fluctv(nbody,ndim),fluct_force(nbody,ndim), meanvel(nphases,ndim)&
         &,varvel(nphases,ndim),meanforce(nphases,ndim)&
         &,meanaccln(nphases,ndim),fluct_accln(nbody,ndim),&
         & acclnstddev(nphases,ndim), fluct_pres_accln(nbody,ndim),&
         & fluct_visc_accln(nbody, ndim), XC_PHYS(nbody,ndim))
    ALLOCATE(AIVJ(nphases,ndim,ndim),VIVJ(nphases,ndim,ndim))

  END SUBROUTINE ALLOCATE_MOVING_PART_RELATED

  SUBROUTINE DEALLOCATE_MOVING_PART_RELATED
    IMPLICIT NONE

    DEALLOCATE(granular_temperature, char_force, meanforcemod,&
         & phase_mass, vivjzioft, vivjetaoft, aivjzioft, aivjetaoft,&
         & source_grantemp, dissip_grantemp,&
         & steady_granular_temp, steady_source_grantemp,&
         & steady_dissip_grantemp,nondim_grantemp&
         &,nondim_source_grantemp,nondim_dissip_grantemp)
    DEALLOCATE(source_pres_grantemp, source_visc_grantemp,&
         & dissip_pres_grantemp, dissip_visc_grantemp, RAD_PHYS)

    DEALLOCATE(fluctv,fluct_force, meanvel,varvel,meanforce,&
         & meanaccln, fluct_accln, acclnstddev, acclnstddev_mod,&
         & fluct_pres_accln, fluct_visc_accln, XC_PHYS)

    DEALLOCATE(AIVJ, VIVJ)

  END SUBROUTINE DEALLOCATE_MOVING_PART_RELATED

  SUBROUTINE open_relevant_files
    IMPLICIT NONE
    CHARACTER(LEN=80) :: FILENAME

    grantempunit = getnewunit(minunitno, maxunitno)
    FILENAME = TRIM(RUN_NAME)//'_granular_temperature'//'.dat'
    OPEN(grantempunit, FILE = FILENAME, STATUS='replace', form&
         &='formatted')
   
    FILENAME = TRIM(RUN_NAME)//'_Aivj'//'.dat'
    aivjunit = getnewunit(minunitno, maxunitno)
    OPEN(aivjunit, FILE = FILENAME, STATUS='replace', form&
         &='formatted')

    FILENAME = TRIM(RUN_NAME)//'_anisotropy'//'.dat'
    anisunit = getnewunit(minunitno, maxunitno)
    OPEN(anisunit, FILE = FILENAME, STATUS='replace', form&
         &='formatted')

    FILENAME = TRIM(RUN_NAME)//'_source_dissip_vs_T'//'.dat'
    source_dissip_unit = getnewunit(minunitno, maxunitno)
    OPEN(source_dissip_unit, FILE = FILENAME, STATUS='replace', form&
         &='formatted')
    !WRITE(*,'(A, 2x, I)')'Unit for source dissipation = .',
    ! source_dissip_unit

    FILENAME = TRIM(RUN_NAME)//'_source_dissip_pres_visc_vs_T'//'.dat'
    source_dissip_pres_visc_unit = getnewunit(minunitno, maxunitno)
    OPEN(source_dissip_pres_visc_unit, FILE = FILENAME, STATUS='replac&
         &e', form='formatted')
    !WRITE(*,'(A, 2x, I)')'Unit for source dissipation pres visc = .',&
    !    & source_dissip_pres_visc_unit
!!$    FILENAME = TRIM(RUN_NAME)//'_eigen_values_vs_t'//'.dat'
!!$    eigenvalue_unit = getnewunit(minunitno, maxunitno)
!!$    OPEN(eigenvalue_unit, FILE = FILENAME, STATUS='replace', form&
!!$         &='formatted')

    FILENAME = TRIM(RUN_NAME)//'_gofr_vs_t'//'.dat'
    gofrunit = getnewunit(minunitno, maxunitno)
    OPEN(gofrunit, FILE = FILENAME, STATUS='replace', form&
         &='formatted')

    FILENAME = TRIM(RUN_NAME)//'_acc_phi_scatter'//'.dat'
    accphiunit = getnewunit(minunitno, maxunitno)
    OPEN(accphiunit, FILE = FILENAME, STATUS='replace', form&
         &='formatted')
  END SUBROUTINE open_relevant_files
  
  SUBROUTINE close_relevant_files
    IMPLICIT NONE

    CLOSE(grantempunit, status = 'keep')
    CLOSE(aivjunit, status = 'keep')
    CLOSE(anisunit, status = 'keep')
    CLOSE(source_dissip_unit, status = 'keep')
    CLOSE(source_dissip_pres_visc_unit, status = 'keep')
!!$    CLOSE(eigenvalue_unit, status = 'keep')
    CLOSE(gofrunit, status = 'keep')
    CLOSE(accphiunit, status = 'keep')
  END SUBROUTINE close_relevant_files

END MODULE post_freely_evolving
