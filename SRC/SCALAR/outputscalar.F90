module outputscalar
#include "../FLO/ibm.h"
  use scalar_data 
  use fftw_interface
  Use dependent_functions
  USE restart_funcs
  USE general_funcs
  USE nlmainarrays, velr=>ubcp, phir=>nlbcp
  USE postproc_funcs
  implicit none 
  Integer :: ix
  !real(prcn), Private,dimension(:,:,:,:), allocatable :: phir, velr
contains
  subroutine output_scal
    implicit none 
    Integer :: i, j, k, count_tmp, runitno, nrbins, ibin, node, strlen
    INTEGER, SAVE :: count_routine=0
    CHARACTER*80 :: SCALANDU_FOR, THETAVSNU3,SCALANDU,SCALONLY,&
         & junk_char, stat, filenameloc 
    LOGICAL :: filexist , isopen
    CHARACTER*10 :: filename2,filename1
    real(prcn), dimension(:), allocatable :: phi_corr, rad_bin
    
    
    nrbins = NINT(real(my,prcn)/1.4d0)
    
    
    !Allocate(phir(mx,my,mz,nspmx), velr(mx,my,mz,ndim))
    count_routine  = count_routine + 1
    IF(count_routine.eq.3) count_routine = 1
    
    if(I_AM_NODE_ZERO) WRITE(*,*) 'COUNT IN SCALAR = ', COUNT_ROUTINE
    !READ(*,*)
    if(I_AM_NODE_ZERO)then
#if PARALLEL
       write (filename1,fmt="('NODE',i2.2,'_',i1)") myid,count_routine
       FILENAMELOC = ""
#else
       WRITE(FILENAME1, '(I1)')count_routine 	
#endif
    end if
    if(I_AM_NODE_ZERO)then
       SCALANDU = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME1)//'.dat&
            &'
       SCALONLY = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME1)//'.dat&
            &'
    else
       SCALANDU = ""
       SCALONLY = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (filename2,fmt="('NODE',i2.2,'_',i1)") node&
               &,count_routine
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME2)//'.dat&
               &'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME2)//'.dat&
               &'

          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(SCALANDU,strlen,node_zero,0,1,decomp_group,status)
       RECV_STRING(SCALONLY,strlen,node_zero,0,1,decomp_group,status)
    end if

        
    INQUIRE(FILE=SCALANDU,EXIST=filexist,OPENED=isopen)
    IF (.NOT.filexist) THEN
       stat = "new"
    ELSEIF(filexist.AND..NOT.isopen) THEN
       stat ="replace"
    ENDIF
    call write_mid_plane (SCALANDU,SCALONLY,stat)
    if(I_AM_NODE_ZERO)then
       SCALANDU_FOR = TRIM(RUN_NAME)//'_SCALANDU_FOR_'&
            &//TRIM(FILENAME1)//'.dat'    
    else
       SCALANDU_FOR = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (filename2,fmt="('NODE',i2.2,'_',i1)") node&
               &,count_routine
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_FOR_'//TRIM(FILENAME2)//'.dat&
               &'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(SCALANDU_FOR,strlen,node_zero,0,1,decomp_group,status)
    end if

    INQUIRE(FILE=SCALANDU_FOR,EXIST=filexist,OPENED=isopen)
    IF (.NOT.filexist) THEN
       stat = "new"
    ELSEIF(filexist.AND..NOT.isopen) THEN
       stat = "replace"
    ENDIF
    call write_forcing_wake(SCALANDU_FOR,stat)
    

!!$    THETAVSNU3 = TRIM(RUN_NAME)//'_THETAVSNU3_'//TRIM(FILENAME1)//'.dat'
!!$    INQUIRE(FILE=THETAVSNU3,EXIST=filexist,OPENED=isopen)
!!$    IF (.NOT.filexist) THEN
!!$       stat = "new"
!!$    ELSEIF(filexist.AND..NOT.isopen) THEN
!!$       stat ="replace"
!!$    ENDIF
    
    !call write_nu_no(THETAVSNU3,stat)

    IF(count_routine.eq.1) count_tmp = 2
    IF(count_routine.eq.2) count_tmp = 1


    runitno = getnewunit(minunitno, maxunitno)
    junk_char = "formatted"
    if(I_AM_NODE_ZERO)then
#if PARALLEL
       write (filename1,fmt="('NODE',i2.2,'_',i1)") myid,count_tmp
       FILENAMELOC = ""
#else
       WRITE(FILENAME1, '(I1)')count_tmp 	
#endif
    end if
    if(I_AM_NODE_ZERO)then
       SCALANDU = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME1)//'.dat&
            &'
       SCALONLY = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME1)//'.dat&
            &'
       SCALANDU_FOR = TRIM(RUN_NAME)//'_SCALANDU_FOR_'//TRIM(FILENAME1)//'.dat&
            &'
    else
       SCALANDU = ""
       SCALONLY = ""
       SCALANDU_FOR = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (filename2,fmt="('NODE',i2.2,'_',i1)") node&
               &,count_tmp
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME2)//'.dat&
               &'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME2)//'.dat&
               &'

          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_FOR__'//TRIM(FILENAME2)//'.dat&
               &'

          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(SCALANDU,strlen,node_zero,0,1,decomp_group,status)
       RECV_STRING(SCALONLY,strlen,node_zero,0,1,decomp_group,status)
       RECV_STRING(SCALANDU_FOR,strlen,node_zero,0,1,decomp_group,status)
    end if
    INQUIRE(FILE=SCALANDU,EXIST=filexist,OPENED=isopen)
    IF(filexist) THEN 
       CALL delete_file(SCALANDU_FOR,junk_char, runitno)
       CALL delete_file(SCALANDU,junk_char, runitno)
       CALL delete_file(SCALONLY,junk_char, runitno)
!!       CALL delete_file(THETAVSNU3,junk_char, runitno)
    end IF

    !ALLOCATE(phi_corr(nrbins), rad_bin(nrbins))
    !CALL scalar_two_point_correlation(phir(1:mx,1:my,1:mz,1:nspmx),phi_corr,rad_bin, nrbins)

    !runitno = getnewunit(minunitno, maxunitno)

    !open(unit=runitno, file = TRIM(RUN_NAME)//'_SCALAR_CORR.dat', form="formatted", status = "unknown")
    !do ibin=1,nrbins
    !   write(runitno,'(4(2x,f12.8))') rad_bin(ibin), rad_bin(ibin)/doml(1), phi_corr(ibin)

    !end do
    !DEALLOCATE(phi_corr,rad_bin)
    !close(runitno,status='keep')


  end subroutine output_scal
  
  subroutine write_mid_plane(scalandu,scalonly,stat)

    Use nlarrays , only : ur1,uf1
    implicit none 
    Integer :: i,j,k, ii, isp
    Real(prcn) :: dist
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: scalandu,scalonly,stat
    LOGICAL :: filexist, isopen
    INTEGER :: count_tmp , scalunit, scalnuunit
    
    
    !    scalunit = getnewunit(minunitno,maxunitno)
!    open(scalunit,FILE=scalonly,form = "formatted",status=stat)
    !    WRITE(scalunit, '("# generated at time = ",g12.5)') t
    scalnuunit = getnewunit(minunitno,maxunitno)
    open(scalnuunit,FILE=scalandu,form = "formatted",status=stat)
    
    WRITE(scalnuunit, '("# generated at time = ",g12.5)') t    
!!$    open(unit=155,file='scalandu.dat',form='formatted',status='unknown')
    write(scalnuunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "UX" ',   &
         &    ' "UY" ',' "UZ" ' ,' "phi" '
    write(scalnuunit,*)'ZONE F=POINT, I=', nx,  ', J=', my
    
!!$    open(unit=156,file='scalonly.dat',form='formatted',status='unknown')
    !write(scalunit,*)'VARIABLES= ',' "X" ',' "Y" ', ' "phi" '
    !write(scalunit,*)'ZONE F=POINT, I=', nx,  ', J=', my
    !Print*,'Hi writing the mid plane file in outputscalar'
    do isp=1,nspmx
       do i=1,nx !mx,1
          
          do j=1,my2,1
             do k=1,mz,1
                uf1(j,k) = phif(i,j,k,isp)
             enddo
          enddo
          
          call ff2cr(uf1,ur1)	 
          do j=1,my,1
             do k=1,mz,1
                phir(i,j,k,isp) = ur1(j,k)+phirmean(isp)
             enddo
          enddo
       enddo
    enddo

    call calc_velreal(u, umean, velr)

    k = mz/2
    do j=1,my
       do i=1,nx !mx
          write(scalnuunit,'(10(3x,g12.5))')real(GLOBAL_INDEX(i)),real(j),velr(i,j,k,1), velr(i,j,k,2),&
               & velr(i,j,k,3), phir(i,j,k,1)!(phir(i,j,k,1)-phistream)/(phisurf&
          ! &-phistream) 
          
       enddo
    enddo

!!$    k = mz/2
!!$    do j=1,my
!!$       do i=1,mx
!!$          !if(j.eq.1.and.i.eq.1) Print*,'foffset = ', foffset
!!$          !if(i.gt.foffset.and.i.lt.(foffset+mxf).and.j.eq.1) then 
!!$          !   Print*, i, (i-foffset-xc(1,1)+radbdy(m))*(one/(two*radbdy(m)))
!!$          !end if
!!$          write(scalunit,'(10(3x,g12.5))') (i-foffset-xc(1,1)+radbdy(1))*(one/(two*radbdy(1))) &
!!$               &,(j-xc(1,2))*(one/(two*radbdy(1))),(phir(i,j,k,1)&
!!$               &-phistream)/(phisurf-phistream)
!!$
!!$
!!$          !write(scalunit,'(10(3x,g12.5))') (i-xc(,real((j-1)),(phir(i,j,k,1)&
!!$          !     &-phistream)/(phisurf-phistream)   
!!$
!!$       enddo
!!$    enddo
    !Print*,'Hi...done writing the mid plane file in outputscalar'
    close(scalnuunit,status='keep')
    close(scalunit,status='keep')
  end subroutine write_mid_plane

  subroutine write_forcing_wake(scalandu_for,stat)
    USE errormesgs
    USE general_funcs
    USE nlmainarrays
    implicit none 
    INTEGER :: i,j,k, unitno
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: scalandu_for,stat
    
    unitno = getnewunit(minunitno, maxunitno)
    OPEN(unitno,FILE=scalandu_for, form="formatted",status=stat)
    
    CALL calc_velreal(u, umean, velr)
    WRITE(unitno, '("# generated at time = ",g12.5)') t
    
!!$    IF (unitno.LT.0) CALL printerror("newunit","ounit")
!!$
!!$    OPEN(unit=unitno,file='scalandu_for.dat',form='formatted',status='re&
!!$         &place')

    write(unitno,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "phi_n" ' ,' &
         &"UX" ',' "UY" ',' "UZ" ' 
    write(unitno,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
    
    if(write_output) then
       do k=1,mz
          do j=1,my
             do i=1,nx
                !print*,'smack my bitch up ', ix
                write(unitno,21) REAL(GLOBAL_INDEX(i)),REAL(j),REAL(k), phir(i,j&
                     &,k,1),velr(i,j,k&
                     &,1) ,velr(i,j,k,2),velr(i,j,k,3) 
             enddo
          enddo
       enddo
    end if
21  FORMAT(10(1xe17.4))
    CLOSE(unitno,status= 'keep')
    
  end subroutine write_forcing_wake
  
  subroutine write_nu_no(THETAVSNU3,stat)
    USE errormesgs
    USE general_funcs
    IMPLICIT NONE 
    INTEGER :: i,j,k,m,l,isp, unitno, iphs
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: thetavsnu3,stat
    LOGICAL :: filexist, isopen

    unitno = getnewunit(minunitno, maxunitno)
    OPEN(unitno,FILE=thetavsnu3, form="formatted",status=stat)
    
!!$    
!!$    IF (unitno.LT.0) CALL printerror("newunit","ounit")
!!$
!!$    OPEN(unit=unitno,file='thetavsNU3.dat',form='formatted',status='re&
!!$         &place')

    DO m=1,nbody
       iphs = part_array(m)%iphs
       nbnd = phase_array(iphs)%nbnd
       nrpr = phase_array(iphs)%nrpr
       
       bndarray => phase_array(iphs)%bndpts
       
       WRITE(unitno,*)'Zone'
       DO isp = 1,nspmx
          DO l=1,nrpr
             IF (bndarray(2,l).GE.zero.AND.bndarray(3,l).EQ.zero) THEN 
                IF(bndarray(1,l).GE.zero) THEN 
                   WRITE(unitno,21)180.-((180.*ATAN(bndarray(2,l)/bndarray(1,l)))/pi)&
                        &,-Nu3(m,l,isp)*two*radbdy(m)*dx
                ELSE
                   WRITE(unitno,21)((180. *ATAN(-bndarray(2,l)/bndarray(1,l)))/pi)&
                        &,-Nu3(m,l,isp)*two*radbdy(m)*dx
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    END DO
21  FORMAT(10(1xe17.4))
    CLOSE(unitno,status= 'keep')

  end subroutine write_nu_no

  SUBROUTINE flow_snapshot
    Use nlarrays , only : ur1,uf1
    USE dem_mod, only : is_mobile, des_pos_new, des_radius
    IMPLICIT NONE
    Integer  :: sunit,i,j,k,m,isp, mark, idim
    INTEGER, SAVE :: zone_count = 0
    LOGICAL, SAVE :: first_time=.TRUE.
    REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
         & mean_vel(ndim), mean_force(ndim), position(ndim)
    CHARaCTER*80 :: FILENAME 
    CHARACTER(LEN=80) :: formfile

    formfile='formatted' 
    
    sunit  = getnewunit(minunitno,maxunitno)
    !if(irestart.eq.1)first_time = .FALSE.
    !CALL calc_velreal(velr)
    if(iscalon.eq.1)then
       do isp=1,nspmx
          do i=1,nx !mx,1
             do j=1,my2,1
                do k=1,mz,1
                   uf1(j,k) = phif(i,j,k,isp)
                enddo
             enddo
             
             call ff2cr(uf1,ur1)	 
             do j=1,my,1
                do k=1,mz,1
                   phir(i,j,k,isp) = ur1(j,k)+phirmean(isp)
                enddo
             enddo
          enddo
       enddo
    else
       do isp=1,nspmx
          do i=1,nx !mx,1
             do j=1,my,1
                do k=1,mz,1
                   phir(i,j,k,isp) = zero
                end do
             end do
          end do
       end do
    end if


    IF(first_time)THEN
       OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted',      &
            &       status='unknown')
       !write(sunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" ',' "UY" ',' "UZ" '!,' "PHI" '
       write(sunit,*)'ZONE T = "', t/t_conv/maxvolfrac, '",'
       write(sunit,*)'DATAPACKING=POINT, I =', nx,  ', J=', my, ', K=', 3    
       first_time = .FALSE.
       k = 1
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz/2
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       FILENAME = TRIM(RUN_NAME)//'_sphr_motion.dat'
       CALL  RUN_TIME_FILE_OPENER(sphrunit,FILENAME, formfile)
    ELSE
       OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted',      &
            &       POSITION='append')
       write(sunit,*)'ZONE T = "', t/t_conv/maxvolfrac, '",'
       write(sunit,*)'DATAPACKING=POINT, I =', nx,  ', J=', my, ', K=', 3
       !write(sunit,*)', VARSHARELIST= ([1-3]=1)'
       k = 1
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz/2
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
    END IF

    close(sunit,status='keep')
!!$    DEALLOCATE(velr)
    
    WRITE(sphrunit,*)'ZONE T= "', t, ' " '
!!$    ucg = DSQRT((rhos/rhof - one)*9.8*dia_phys)
    do idim = 1, ndim
       mean_force(idim) = SUM(force(1:nbody,idim))/real(nbody,prcn)
       mean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
    end do

    DO m=1,nbody
!!$       WRITE(sphrunit,'(4(2x,f12.8))')  t*ucg/dia_phys,velbdy(m,1:ndim)/ucg
!!$       WRITE(*,'(4(2x,f12.8))')  t*ucg/dia_phys,velbdy(m,1:ndim)/ucg

       fluct_vel(1:ndim) = velbdy(m,1:ndim)-mean_vel(1:ndim)
       fluct_force(1:ndim) = force(m,1:ndim)-mean_force(1:ndim)
       mark = -1
       if(DOT_PRODUCT(fluct_force(1:ndim),fluct_vel(1:ndim)).gt.zero) mark = 1

       do idim = 1, ndim
          position(idim) = XC(m,idim)+frame_pos(idim)
          if(position(idim).lt.one) position(idim) = position(idim)+&
               & real(my,prcn)
          if(position(idim).ge.real(my+1,prcn))position(idim) =&
               & position(idim)- real(my,prcn)
       end do

       WRITE(sphrunit,'(10(2x,f12.8))')  position(1), position(2), position(3), radbdy(m),radbdy(m),real(mark)
    enddo
!!$    close(sphrunit,status='keep')
    !Write(*,*) 'WRITITNG SNAPSHOT FILE, ZONE COUNT = ', zone_count
  END SUBROUTINE flow_snapshot

#if 0
  SUBROUTINE flow_snapshot2
    Use nlarrays , only : ur1, uf1
	 use nlmainarrays, only : ubcp
    USE dem_mod, only : is_mobile, des_pos_new, des_radius
    IMPLICIT NONE
    Integer  :: sunit,i,j,k,l,m,isp, mark, idim
    INTEGER, SAVE :: zone_count = 0
    LOGICAL, SAVE :: first_time=.TRUE.
    REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
         & mean_vel(ndim), mean_force(ndim), position(ndim)
    CHARaCTER*80 :: FILENAME1, filename2, filename3
    integer, save :: sphrunit1, sphrunit2, sphrunit3
    CHARACTER(LEN=80) :: formfile

	real(8), allocatable :: out_arr(:,:,:,:), trans_buf(:)
	integer :: node_num, iproc
	integer :: jj, j1, j2, j3
	real(8) :: tmp, tmp1, tmp2, tmp3
	integer :: iphs, part_start, part_end
	logical :: filexist
	

	formfile='formatted' 

	j1=1
	j2=my/2
	j3=my

	if (I_AM_NODE_ZERO) write (*,*) "IN FLOW_SNAPSHOT"

#if PARALLEL
	if (I_AM_NODE_ZERO) then
		allocate(out_arr(mx1,3,mz,ndim+1))

		do jj=1, 3
			if (jj==1) then
				j=j1
			elseif (jj==2) then
				j = j2
			elseif (jj==3) then
				j = j3
			endif

			! velocity fluctuations for node zero
			do idim=1, ndim
				do i=1, nx
					do k=1, mz
						out_arr(i,jj,k,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
!							if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
					enddo
				enddo
			enddo
		enddo

		do iproc=1,nproc-1
			node_num = mz*(ends(iproc)-starts(iproc)+1)*3
			allocate(trans_buf(node_num))

			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				! collecting velocity fluctuations from other processes
				call mpi_recv(trans_buf(1),node_num,mpi_double_precision,iproc,iproc,decomp_group,status,err_code)

				l=0
				do idim=1, ndim
					do k=1, mz
						do i=starts(iproc),ends(iproc)
							l=l+1
							out_arr(i,jj,k,idim) = trans_buf(l)
						enddo
					enddo
				enddo
			enddo
			deallocate(trans_buf)
		enddo
	else
		! recieving velocity fluctuations from node zero
		node_num=mz*nx*3
		allocate(trans_buf(node_num))

		do jj=1, 3
			if (jj==1) then
				j=j1
			elseif (jj==2) then
				j = j2
			elseif (jj==3) then
				j = j3
			endif

			l=0
			do idim=1,ndim
				do k=1, mz
					do i=1, nx
						l=l+1
						trans_buf(l) = ubcp(i,j,k,idim) !-ufmean(idim)
!						if (.not.fluid_atijk(i,j,k)) trans_buf(l) = 0d0
					enddo
				enddo
			enddo

			call mpi_send(trans_buf(1),node_num,mpi_double_precision,node_zero,myid,decomp_group,err_code)
		enddo
		deallocate(trans_buf)
	endif
#else
	allocate(out_arr(mx1,3,mz,ndim+1))
	do idim=1, ndim
		do k=1, mz
			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				do i=1, mx1
					out_arr(i,jj,k,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
!					if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
				enddo
			enddo
		enddo
	enddo
#endif

	if (I_AM_NODE_ZERO) then
		write (*,*) "GENERATING THE SNAPSHOT OF THE FIELD"

		do k=1, mz
			do j=1, 3
				do i=1, mx1
					tmp = 0d0
					do idim=1, ndim
						tmp = tmp + (out_arr(i,j,k,idim)-ufmean(idim)) * (out_arr(i,j,k,idim)-ufmean(idim))
					enddo
					out_arr(i,j,k,4) = tmp / umeanslip**2
					out_arr(i,j,k,1:ndim) = out_arr(i,j,k,1:ndim) / umeanslip
				enddo
			enddo
		enddo

      sunit     = 30
		sphrunit1 = 31 !getnewunit(minunitno,maxunitno)
		sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
		sphrunit3 = 33 !getnewunit(minunitno,maxunitno)

		FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'
		filename2 = TRIM(RUN_NAME)//'_sphr_motion_act1.dat'
		filename3 = TRIM(RUN_NAME)//'_sphr_motion_act2.dat'


		inquire (file=trim(filename1), exist=filexist)
		if (.not.filexist) then
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='replace')

!			sphrunit  = 31 !getnewunit(minunitno,maxunitno)
!			sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
!			sphrunit3 = 33 !getnewunit(minunitno,maxunitno)
!
!			FILENAME = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'
!			filename2 = TRIM(RUN_NAME)//'_sphr_motion_act1.dat'
!			filename3 = TRIM(RUN_NAME)//'_sphr_motion_act2.dat'

!			CALL  RUN_TIME_FILE_OPENER(sphrunit,FILENAME, formfile)
!			CALL  RUN_TIME_FILE_OPENER(sphrunit2,FILENAME2, formfile)
!			CALL  RUN_TIME_FILE_OPENER(sphrunit3,FILENAME3, formfile)

			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='unknown')
!			OPEN(unit = sphrunit2,file=TRIM(FILENAME2),form='formatted', status='unknown')
!			OPEN(unit = sphrunit3,file=TRIM(FILENAME3),form='formatted', status='unknown')
		ELSE
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='old', position='append')

			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='old', position='append')
!			OPEN(unit = sphrunit2,file=TRIM(FILENAME2),form='formatted', status='old', position='append')
!			OPEN(unit = sphrunit3,file=TRIM(FILENAME3),form='formatted', status='old', position='append')
		endif

		write(sunit,*)'ZONE T = "', t/t_conv/maxvolfrac, '",'
		write(sunit,*)'DATAPACKING=POINT, I =', mx1,  ', J=', 3, ', K=', mz


		do k=1, mz
			do jj=1, 3
				if (jj==1) then
					j=j1
				elseif (jj==2) then
					j = j2
				elseif (jj==3) then
					j = j3
				endif

				do i=1, mx1
					write(sunit,"(3i6,1d15.7)") i, j, k, out_arr(i,jj,k,4) !out_arr(i,jj,k,1), out_arr(i,jj,k,2), out_arr(i,jj,k,3), out_arr(i,jj,k,4)
				enddo
			enddo
		enddo
		close(sunit,status='keep')

		WRITE(sphrunit1,*)'ZONE T= "', t/t_conv/maxvolfrac, ' " '

		do iphs=1, nphases
			part_start = phase_array(iphs)%pstart
			part_end   = phase_array(iphs)%pend

			do m=part_start, part_end
				write (sphrunit1,"(6d15.7, 1i6)")  xc(m,:), radbdy(m), sqrt(dot_product(velbdy(m,:),velbdy(m,:)))/umeanslip, &
					& sqrt(dot_product(phase_array(iphs)%mean_spec_vel(:),phase_array(iphs)%mean_spec_vel(:)))/umeanslip, iphs
			enddo
		enddo

!		do m=1, nbody
!			if (m==1) then
!				write (sphrunit1,"(5d15.7)")  xc(m,:), radbdy(m), one
!			else
!				write (sphrunit1,"(5d15.7)")  xc(m,:), radbdy(m), zero
!			endif
!			write (sphrunit1,"(5d15.7)")  xc(m,:), radbdy(m), sqrt(dot_product(velbdy(m,:),velbdy(m,:))) !/umeanslip
!		enddo




!		WRITE(sphrunit2,*)'ZONE T= "', t, ' " '
!		WRITE(sphrunit3,*)'ZONE T= "', t, ' " '


		close(sphrunit1)
!		close(sphrunit2)
!		close(sphrunit3)
!		stop

!		do idim = 1, ndim
!			mean_force(idim) = SUM(force(1:nbody,idim))/real(nbody,prcn)
!			mean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
!		end do

!		DO m=1,nbody
!			fluct_vel(1:ndim) = velbdy(m,1:ndim)-mean_vel(1:ndim)
!			fluct_force(1:ndim) = force(m,1:ndim)-mean_force(1:ndim)
!			mark = -1
!			if(DOT_PRODUCT(fluct_force(1:ndim),fluct_vel(1:ndim)).gt.zero) mark = 1
!			do idim = 1, ndim
!				position(idim) = XC(m,idim) !+frame_pos(idim)
!				if(position(idim).lt.one) position(idim) = position(idim)+ real(my,prcn)
!				if(position(idim).ge.real(my+1,prcn)) position(idim) = position(idim)- real(my,prcn)
!			end do

!			tmp1 = j1*dy
!			tmp2 = j2*dy
!			tmp3 = j3*dy
!
!			if (abs(xc(m,2)-j1)<=radbdy(m).or.abs(xc(m,2)-j3)<=radbdy(m)) then
!				WRITE(sphrunit2,"(4d15.7)")  xc(m,:), radbdy(m) !,radbdy(m),real(mark)
!			elseif (abs(xc(m,2)-j2)<=radbdy(m)) then
!				WRITE(sphrunit3,"(4d15.7)")  xc(m,:), radbdy(m) !,radbdy(m),real(mark)
!			else
!				WRITE(sphrunit1,"(4d15.7)")  xc(m,:), radbdy(m) !,radbdy(m),real(mark)
!			endif
!		enddo
!		close (sphrunit1)
!		close (sphrunit2)
!		close (sphrunit3)

		deallocate(out_arr)
	endif
  END SUBROUTINE flow_snapshot2









  SUBROUTINE flow_snapshot3
    Use nlarrays , only : ur1, uf1
	 use nlmainarrays, only : ubcp
    USE dem_mod, only : is_mobile, des_pos_new, des_radius
    IMPLICIT NONE
    Integer  :: sunit,i,j,k,l,m,isp, mark, idim
    INTEGER, SAVE :: zone_count = 0
    LOGICAL, SAVE :: first_time=.TRUE.
    REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
         & mean_vel(ndim), mean_force(ndim), position(ndim)
    CHARaCTER*80 :: FILENAME1, filename2, filename3
    integer, save :: sphrunit1, sphrunit2, sphrunit3
    CHARACTER(LEN=80) :: formfile

	real(8), allocatable :: out_arr(:,:,:,:), trans_buf(:)
	integer :: node_num, iproc
	integer :: kk, k1, k2, k3
	real(8) :: tmp, tmp1, tmp2, tmp3
	integer :: iphs, part_start, part_end
	logical :: filexist
	

	formfile='formatted' 

	k1=1
	k2=mz/2
	k3=mz

	if (I_AM_NODE_ZERO) write (*,*) "IN FLOW_SNAPSHOT"

#if PARALLEL
	if (I_AM_NODE_ZERO) then
		allocate(out_arr(mx1,my,3,ndim+1))

		do kk=1, 3
			if (kk==1) then
				k=k1
			elseif (kk==2) then
				k = k2
			elseif (kk==3) then
				k = k3
			endif

			! velocity fluctuations for node zero
			do idim=1, ndim
				do i=1, nx
					do j=1, my
						out_arr(i,j,kk,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
!							if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
					enddo
				enddo
			enddo
		enddo

		do iproc=1,nproc-1
			node_num = my*(ends(iproc)-starts(iproc)+1)*3
			allocate(trans_buf(node_num))

			do kk=1, 3
				if (kk==1) then
					k = k1
				elseif (kk==2) then
					k = k2
				elseif (kk==3) then
					k = k3
				endif

				! collecting velocity fluctuations from other processes
				call mpi_recv(trans_buf(1), node_num, mpi_double_precision, iproc, iproc, decomp_group, status, err_code)

				l=0
				do idim=1, ndim
					do j=1, my
						do i=starts(iproc),ends(iproc)
							l=l+1
							out_arr(i,j,kk,idim) = trans_buf(l)
						enddo
					enddo
				enddo
			enddo
			deallocate(trans_buf)
		enddo
	else
		! recieving velocity fluctuations from node zero
		node_num=my*nx*3
		allocate(trans_buf(node_num))

		do kk=1, 3
			if (kk==1) then
				k=k1
			elseif (kk==2) then
				k = k2
			elseif (kk==3) then
				k = k3
			endif

			l=0
			do idim=1, ndim
				do j=1, my
					do i=1, nx
						l=l+1
						trans_buf(l) = ubcp(i,j,k,idim) !-ufmean(idim)
!						if (.not.fluid_atijk(i,j,k)) trans_buf(l) = 0d0
					enddo
				enddo
			enddo

			call mpi_send(trans_buf(1),node_num,mpi_double_precision,node_zero,myid,decomp_group,err_code)
		enddo
		deallocate(trans_buf)
	endif
#else
	allocate(out_arr(mx1,my,3,ndim+1))
	do idim=1, ndim
		do j=1, my
			do kk=1, 3
				if (kk==1) then
					k=k1
				elseif (kk==2) then
					k = k2
				elseif (kk==3) then
					k = k3
				endif

				do i=1, mx1
					out_arr(i,j,kk,idim) = ubcp(i,j,k,idim) !-ufmean(idim)
!					if (.not.fluid_atijk(i,j,k)) out_arr(i,jj,k,idim) = 0d0
				enddo
			enddo
		enddo
	enddo
#endif

	if (I_AM_NODE_ZERO) then
		write (*,*) "GENERATING THE SNAPSHOT OF THE FIELD"

		do k=1, 3
			do j=1, my
				do i=1, mx1
					tmp = 0d0
					do idim=1, ndim
						tmp = tmp + (out_arr(i,j,k,idim)-ufmean(idim)) * (out_arr(i,j,k,idim)-ufmean(idim))
					enddo
					out_arr(i,j,k,ndim+1) = tmp / umeanslip**2
					out_arr(i,j,k,1:ndim) = out_arr(i,j,k,1:ndim) / umeanslip
				enddo
			enddo
		enddo

      sunit     = 30
		sphrunit1 = 31 !getnewunit(minunitno,maxunitno)
		sphrunit2 = 32 !getnewunit(minunitno,maxunitno)
		sphrunit3 = 33 !getnewunit(minunitno,maxunitno)

		FILENAME1 = TRIM(RUN_NAME)//'_sphr_motion_pas.dat'

		inquire (file=trim(filename1), exist=filexist)
		if (.not.filexist) then
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='replace')
			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='unknown')
		ELSE
			OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted', status='old', position='append')
			OPEN(unit = sphrunit1,file=TRIM(FILENAME1),form='formatted', status='old', position='append')
		endif

		write(sunit,*)'ZONE T = "', t/t_conv, '",'
		write(sunit,*)'DATAPACKING=POINT, I =', mx1,  ', J=', my, ', K=', 3

		do kk=1, 3
			if (kk==1) then
				k=k1
			elseif (kk==2) then
				k = k2
			elseif (kk==3) then
				k = k3
			endif

			do j=1, my
				do i=1, mx1
					write(sunit,"(3i6,1d15.7)") i, j, k, out_arr(i,j,kk,4) !out_arr(i,jj,k,1), out_arr(i,jj,k,2), out_arr(i,jj,k,3), out_arr(i,jj,k,4)
				enddo
			enddo
		enddo
		close(sunit,status='keep')

		WRITE(sphrunit1,*)'ZONE T= "', t/t_conv, ' " '

		do iphs=1, nphases
			part_start = phase_array(iphs)%pstart
			part_end   = phase_array(iphs)%pend

			do m=part_start, part_end
				write (sphrunit1,"(5d15.7)")  xc(m,:), radbdy(m), sqrt(dot_product(velbdy(m,:),velbdy(m,:)))/umeanslip
			enddo
		enddo
		close(sphrunit1)
		deallocate(out_arr)
	endif
  END SUBROUTINE flow_snapshot3
#endif
end module outputscalar
