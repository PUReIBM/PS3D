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

MODULE scalar_init
#include "../FLO/ibm.h"
  USE scalar_data
  USE fftw_interface
  USE nlarrays, ONLY  : phif1=>uf1, phir1=>ur1
  !  Use test_flux 
  IMPLICIT NONE 
  
  REAL(prcn) ::  xl,y,z,temp
  
  INTEGER ::  diff2
CONTAINS 
  SUBROUTINE initscal
    IMPLICIT NONE   
    INTEGER :: i,j,k,isp
    REAL*8 :: lchar 

    nx = xend-xstart+1
    if(I_AM_NODE_ZERO) write(*,'(A)')'INITIALIZING AND ALLOCATING THE SCALAR DATA '
    
    !call alloc_scalar_mem
    !Initialize the surface temperature and stream temp.
    
    IF(NSPMX.GT.3) THEN
       if(I_AM_NODE_ZERO)then
          PRINT*,'USING THE FLOW FR ARRAY, WHICH IS SIZE 3'
          PRINT*,'MAKE SEPARATE ARRAY IN SCALAR_BCSET'
          WRITE(OUNIT,*)'USING THE FLOW FR ARRAY, WHICH IS SIZE 3'
          WRITE(OUNIT,*)'MAKE SEPARATE ARRAY IN SCALAR_BCSET'
       end if
       PARALLEL_FINISH()
       STOP
    end IF
    
    !inner_itn = .False.
    
    !gamma(1) = kf/(rhof*cpf)
    
    gamma(1) = vis/Pr_or_Sc_nu
    
    !Lchar = (pi/(6.d0*maxvolfrac))**(one/three)
    lchar = char_length !dchar!*(one-maxvolfrac)

    
    t_diff = (lchar*lchar)/MAXVAL(gamma(1:nspmx))
    if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.6)') 'Pr = tdiff/tvis = ', t_diff/t_vis
    
    if(t_diff-t_min.LT.small_number) then 
       t_min = t_diff
       diff_ts = .true.
       
       lchar = ((one-maxvolfrac)*doml(1)*doml(2)*doml(3))**(one/three)
       
       if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.6)') 'LCHAR FOR VIS EXT LEN = ', LCHAR
       tendused = 0.2*tend*lchar*lchar/MAXVAL(gamma(1:nspmx))
       
       
       !dt_tmp_diff = (vfl*cfl*lchar*lchar*(one-maxvolfrac))/MAXVAL(gamma(1:nspmx))
    end if
    
    lchar = dx 
    dt_tmp_diff = (lchar*lchar*(one-maxvolfrac))/MAXVAL(gamma(1:nspmx))
    dt_tmp = dt_tmp_diff
    if(I_AM_NODE_ZERO)  WRITE(*,'(A,3(1x,g12.5))') 'T_CONV, T_VIS, T_DIFF = ',t_conv,t_vis, T_DIFF
    
    !PRINT*,'SIZE = ', SIZE(phisurfall,1), SIZE(phisurfall,2), nbody
    phisurfall(1:nbody,1) = phisurf
    
    
    !dt_tmp_diff = (vfl*cfl*lchar*lchar)/MAXVAL(gamma(1:nspmx))
    
    IF(iscal_restart.eq.0) then 
       
       sourcesink(:) = 0.0d0
       fphirmean = zero
       
       nu_old = zero
       
       nu_error_array(:) = one
       nu_error_hist = one
       nu_error = one 
    end IF
    

    if(I_AM_NODE_ZERO) write(*,*) 'GAMMA HAS BEEN SET SO AS TO MAKE Pr =',  Pr_or_Sc_nu
    if(I_AM_NODE_ZERO) write(ounit,*) 'GAMMA HAS BEEN SET SO AS TO MAKE Pr',  Pr_or_Sc_nu
    if(I_AM_NODE_ZERO) write(*,'(2x,A6,g12.5,/,2x,A6,i3,/,2x,A6,g12.5,/,2x,A6,g12.5,/,&
         & 2x,A6,g12.5,/,2x,A10,g12.5,/)') &
         &'gamma=',gamma(1), &
         &'nsp=',nspmx, &
         &'rhos =',rhos, &
         &'cps =',cps, &
         &'ks = ', ks,  &
         & 'Prandtl # = ', Pr_or_Sc_nu
    !u(:,:,:,:)= czero!only for the testing 
    
    !----------------------------------------------------------------------
    diff = 10
    diff2 = 10
    nxi = (30-diff)
    nxf = (30+diff)
    nyi = (my/2-diff2)
    nyf = (my/2+diff2)
    nzi = (mz/2-diff2)
    nzf = (mz/2+diff2)
    xis = (nxi-1)*dx
    xf = (nxf-1)*dx
    yi = (nyi-1)*dy
    yf = (nyf-1)*dy
    zi = (nzi-1)*dz
    zf = (nzf-1)*dz
    !print*,'xf=',xf,'xis=',xis,'diff=',diff
    cons=1.d0
    uo=1.d0
    !       Load the initial Scalar field
    !WRITE(*,*)'Reading initial scalar field'
    
    IF(iscal_restart.EQ.0)THEN
       if(setphimean) then
          phirmean(1:nspmx) = (one-maxvolfrac)*phistream+(maxvolfrac)*phisurf
          
       else 
          !phirmean(1:nspmx) = phistream
          phirmean(1:nspmx) = (one-maxvolfrac)*phistream+(maxvolfrac)*phisurf
          phimean_des = (one-maxvolfrac)*phistream+(maxvolfrac)*phisurf
       end if
       
       !phirmean(1:nspmx) = phistream
       !WRITE(*,*)'Begin..transforming real velocity to Fourier space'
!!$       DO isp=1,nspmx
!!$          DO i=1,mx,1
!!$             DO k=1,mz,1
!!$                DO j=1,my,1
!!$                   phir1(j,k) = zero
!!$                   
!!$                END DO
!!$             END DO
!!$             CALL ff2rc(phir1(:,:),phif(i,:,:,isp))
!!$          END DO
!!$       END DO
       !ofr(:,:,:,:) = zero 
    END IF
    
    phimean_des = (one-vol_frac1)*phistream+(vol_frac1)*phisurf

    !Initialize boundary conditions for scalars in Fourier Space
    !phif(1,1,1,1) = dcmplx(4.0,zero)
    
!!$    DO isp = 1, nspmx,1
!!$       DO k=1,mz
!!$          DO j=1,my2
!!$             phiin(j,k,isp) = phif(1,j,k,isp)
!!$             !		 phiin(j,k,isp)=czero
!!$             phiout(j,k,isp) = phif(mx,j,k,isp) 
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
    IF (sourcepresent) THEN 
       !READ(*,*)
#if PARALLEL
       if(I_AM_NODE_ZERO) PRINT*,'ADDING SOURCE NOT AVAILABLE IN PARAL&
            &LEL VERSION. PLEASE PARALLELIZE. EXITING'
       PARALLEL_FINISH()
       STOP
#else
       CALL add_source
#endif
    end IF

    !gamma = zero 
  END SUBROUTINE initscal

  SUBROUTINE alloc_scalar_mem
  USE global_data
  USE scalar_data
#if PARALLEL
  USE nlarrays, Only : uatminus1, uatnxp2
#endif
  IMPLICIT NONE
  INTEGER :: i,j,k
  
  !ALLOCATE COMPLEX ARRAYS 
  !PRINT*,'IN ALLOC SCALAR MEM: NBODY = ', NBODY

#if !PARALLEL
  ALLOCATE(phif(nx+1,my2,mz,nspmx))
  ALLOCATE(nlphif(nx+1,my2,mz,nspmx),onlphif(nx+1,my2,mz,nspmx))
#else
  ALLOCATE(phif(0:nx+1,my2,mz,nspmx))
  ALLOCATE(nlphif(nx,my2,mz,nspmx),onlphif(nx,my2,mz,nspmx))
  ALLOCATE(uatminus1(my,mz), uatnxp2(my,mz))
#endif

  
  !ALLOCATE(gradphi(mx,my,mz,3,nspmx))
  ALLOCATE(sourcesink(nspmx))
  ALLOCATE( flux_global(nspmx))
  ALLOCATE(flux_global2(nspmx))
  !ALLOCATE(surf_scal_value(nbnd,nrpr,nspmx), Nu3(nbnd,nrpr,nspmx))
  sourcesink = zero 
  do k = 1, mz
     do j = 1, my2
#if PARALLEL
        phif(0,j,k,:) = czero
#endif
        
        do i = 1,nx+1
#if PARALLEL
           if(i.le.nx)then
#endif
              phif(i,j,k,:) = czero
              !ffphi(i,j,k,:) = czero
              nlphif(i,j,k,:) = czero
              onlphif(i,j,k,:) = czero
#if PARALLEL
           end if
#endif
           !gradphi(i,j,k,:,:) = zero 
        end do
#if PARALLEL
        phif(nx+1,j,k,:) = czero
#endif
     end do
  end do
    
!!$  phiin = czero 
!!$  phiout = czero
  allocate(gamma(nspmx))

  
  ALLOCATE(nu_error_array(nerr_steps))
  ALLOCATE(phi_fluid_mean(nspmx))
  ALLOCATE(phi_solid_mean(nspmx))
  ALLOCATE(phimodmean(nspmx))
  
  end SUBROUTINE alloc_scalar_mem
  
  
  SUBROUTINE add_source 
    IMPLICIT NONE 
    INTEGER :: i,j,k,isp
    
    ALLOCATE(sourcephif(mx,my2,mz,nspmx))
    ALLOCATE(sourcephi(mx,my,mz,nspmx))
    OPEN(unit=33,file='source.dat',status='unknown')
    !Initialize a source term for the scalar field
    
    WRITE(33,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "SOU" ',' "IC" '
    WRITE(33,*)'ZONE F=POINT, I=',mx , ', J=',my, ', K=',mz
    DO isp = 1,nspmx
       
       !          do k = nzi,nzf
       DO k = 1,mz
          !	      do j=nyi,nyf
          DO j = 1,my
             !                 do i = nxi,nxf
             DO i = 1,mx
                xl = (i-1)*dx
                y=(j-1)*dy
                z=(k-1)*dz
                !sourcephi(i,j,k,isp) = phistream!sin((2*pi/doml(1))*xl)
                sourcephi(i,j,k,isp) = sin((2.d0*pi/doml(2))*y)+sin((2*pi/doml(1))*xl)
                !sourcephi(i,j,k,isp)=cons*sin(pi*(xl-xis)/(xf
                !PRINT*, 'SORUCE + ',                 sourcephi(i,j,k,isp)
                !-xis))*sin(pi*(y-yi)/(yf-yi))*sin(pi*(z-zi)/(zf
                !-zi)) 
                !sourcephi(i,j,k,isp)=cons*cos(pi*(y-yi)/(yf-yi))
                !sourcephi(i,j,k,isp)=cons*cos(pi*(z-zi)/(zf-zi))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !	do k = 1,mz
    !	   do j = 1,my
    !	      do i = 1,mx
    !		 write(33,*) (i-1)*dx,(j-1)*dy,(k-1)*dz,
    !	1	      sourcephi(i,j,k,1),
    !	1	      phireal(i,j,k,1)
    !		 print*,sourcephi(i,j,k,1)
    !	      enddo
    !	   enddo
    !	enddo
    DO isp = 1,nspmx
       DO i=1,mx
          DO k=1,mz,1
             DO j=1,my,1
                phir1(j,k)=sourcephi(i,j,k,isp)

             ENDDO
          ENDDO
          
          CALL ff2rc(phir1,phif1)
          DO k=1,mz
             DO j=1,my2
                sourcephif(i,j,k,isp)=phif1(j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    PRINT*,'source AVG = ',SUM(sourcephif(1:mx1,1,1,1))/(mx1)
    
    !DEALLOCATE(sourcephi)
  END SUBROUTINE add_source

SUBROUTINE init_params_scal
  IMPLICIT NONE 
  CALL set_default_scal
  
  CALL read_nmls_scal 
  
  !call write_nmls_scal
  
  CALL write_input_scal(ounit)
END SUBROUTINE init_params_scal

SUBROUTINE set_default_scal
  
  IMPLICIT NONE
  INTEGER :: isp 
  nspmx=1
  phistream = 1.d0
  phisurf = 0.d0
  DO isp = 1, nspmx
     gamma(isp) =  2.8169E-03
  END DO
  sourcepresent = .FALSE.
  zero_flow_in_solid = .TRUE.
END SUBROUTINE set_default_scal


SUBROUTINE read_nmls_scal
  
  Use errormesgs
  Use general_funcs
  implicit none 
  LOGICAL:: nmlexist
  INTEGER:: unitno, ierr
  
   unitno = getnewunit(minunitno,maxunitno)
  OPEN(unit=unitno,file="scalparam.in",form="formatted",DELIM='APOSTROP&
       &HE', IOSTAT=ierr)
  PRINT*,'READING THE SCALAR NAMELIST FILE scalparam.in'
  !nmlexist = positionnml(unitno,"Properties")
  READ(unitno,NML=scal_propt)
  
  close(unitno,STATUS='keep')

END SUBROUTINE read_nmls_scal

SUBROUTINE write_nmls_scal !Just a debug to see how the namelists are
  ! written out 
  implicit none 
  LOGICAL:: nmlexist
  INTEGER:: unitno

  unitno = 101     ! have to use arbitrary unit number
  OPEN(unit=unitno,file="scalar.opt",form="formatted",DELIM='APOSTROPHE')
  !PRINT*,'HI... writing out the namelists file'
  !nmlexist = positionnml(unitno,"Properties")
  
  WRITE(unitno,NML=scal_propt)
  CLOSE(unitno)
  
END SUBROUTINE write_nmls_scal
  
SUBROUTINE write_input_scal(fileno)
  
  Use errormesgs
  Use general_funcs
  implicit none 
  Integer, Intent(in):: fileno
  !Print*,'Hi ... writing namelists to outputfile...filno=',fileno!DEBUG
  call separator(fileno,43,'*')
  write(fileno,*)'INPUT PARAMETERS FOR THE SCALAR'
  call separator(fileno,43,'*')
  call blankline(fileno)
  write(fileno,NML=scal_propt)
end SUBROUTINE write_input_scal
 

SUBROUTINE save_unformatted_scal
  implicit none 
  
  WRITE(*,*) '-------------------------------------'
  WRITE(*,*) 'WRITING SCALAR UNFORMATTED RESTART DATA'
  if(iscalon.eq.1) then
     Open(111, file='scal.rst', form='unformatted')
     write(111)phif
     write(111)nlphif
     write(111)onlphif
     write(111) phisurfall(1:nbody,1:nspmx)
     write(111) phirmean(1:nspmx), fphirmean(1:nspmx)
     close(111)
  endif
  
  WRITE(*,*) '-------------------------------------'
end SUBROUTINE save_unformatted_scal

SUBROUTINE read_unformatted_scal
  Implicit none
  WRITE(*,*) '-------------------------------------'
  WRITE(*,*) 'IN SCALAR READ UNFORMATTED RESTART'
  
  if(iscalon.eq.1) then
     Open(111, file='scal.rst', form='unformatted')
     read(111)phif
     read(111)nlphif
     read(111)onlphif
     read(111) phisurfall(1:nbody,1:nspmx)
     read(111)phirmean(1:nspmx),  fphirmean(1:nspmx)
     close(111)
    endif

    PRINT*,'phirmean and fphirmean = ', phirmean, fphirmean 
    WRITE(*,*) '-------------------------------------'
  
  END SUBROUTINE read_unformatted_scal

END MODULE scalar_init
    
