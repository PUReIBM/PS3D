MODULE spectra_3d
  USE precision  
  USE constants  
  USE general_funcs
  USE global_data
  USE scalar_data
  USE dependent_functions
  USE post_global_data
  USE nlmainarrays, ur_tot=>ubc, arr=>nlbc, misvel=>onlbc
  use fftw_interface

  implicit none

  PRIVATE

  TYPE GRID
     REAL(prcn) :: dx, dy, dz
     INTEGER :: cx, cy, cz ! Number of grid nodes
     REAL(prcn), DIMENSION(:), POINTER :: XE, YN, ZT
     REAL(prcn) :: Ener, Ener_res, Ener_cross
     !REAL(prcn), DIMENSION(:,:,:), POINTER :: meannum, meannumsq
     !REAL(prcn), DIMENSION(:,:,:), POINTER :: meanvol, meanvolsq
  END TYPE GRID
  TYPE(GRID) :: coarse_grid
  REAL(prcn) :: meanuf(ndim), misavg_mean_vel(ndim), parallel_dir(ndim)
  
  PUBLIC :: calc_residual_energy, calc_fluid_scalar_stats,&
       & calc_filtered_mis_avg_fields, calc_mis_avg_vel_field

CONTAINS 
  SUBROUTINE calc_spec_3d
    implicit none
    double complex, dimension(my2,my,mz,ndim) :: turb_uf
    double complex :: wx(my2), wy(my), wz(mz)
    double precision, allocatable, dimension(:) :: E_spectra, karray
    !integer :: particle_pos(0:mx+1,my,mz)
    double precision :: k0, kmax, tmp, kmode, lx, ly, lz, tmp1, lenscale, tketotal, kmin 
    double precision :: umean_int(ndim), tke(ndim),  tke_r(ndim)
    double precision :: kmaxeta, eta, dkappa
    integer :: i,j,k,idim, count,turb_lx, ibody, istat, count_real, count_comp,kindex, kindex_max

    integer :: iloc, jloc, kloc, xlo, xhi, ylo, yhi, zlo, zhi, inpos, outpos, nkbins
    double precision :: dist, energy_bal, epsf, ktot, kmode_max,&
         & dxlength, avg_phii, decimal, E_spectrum_function
    REAL(prcn) :: mean_energy

    ly = doml(2)
    lz = doml(3)
    lx = doml(1)
    print *,'lx,ly,lz',lx,ly,lz

    !mean_energy = (ucharmod/(1-vol_frac1))**2.d0
    
!!$    do i=1,my
!!$       do idim=1,ndim
!!$          call ff2cr(u(i,:,:,idim), ur_tot(i,:,:,idim))
!!$       end do
!!$    end do


    do i=1, my/2+1
       wx(i)=dcmplx(zero,2.0*pi*dble(i-1)/dble(lx))
    enddo
    
    do i=1,my/2
       wy(i) = dcmplx(zero,2*pi*dble(i-1)/dble(ly))
       wy(my+1-i) = dcmplx(zero,-1.d0*2*pi*dble(i)/dble(ly))
    enddo
    
    do i=1,mz/2
       wz(i)=dcmplx(zero,2.0*pi*dble(i-1)/dble(lz))
       wz(mz+1-i)=dcmplx(zero,-1.d0*2*pi*dble(i)/dble(lz))
    enddo
    
    
    kmaxeta = 1.0
    k0 = aimag(wx(2))
    !kmax = sqrt(2.0)/3.d0*my*k0
    kmax = my*k0/two ! Smallest wave that can be resolved is of
    ! wavelength size 2*dx (not dx). 
    eta = kmaxeta/kmax
    !dkappa shud be in the radial/spherical coordinate system. THe above wavenumbers are defined
    !on a cartesian grid and dkappa is chosen in the radial system
    dkappa = dsqrt(3.d0)*two*pi/lx
    Write(*,*)'dkappa = ', dkappa
    
    kmin = two*pi/lx
    Write(*,*)'Minimum wavenumber = ', kmin
 
    nkbins = INT((kmax-kmin)/dkappa)
    !nkbins = INT(kmax/dkappa)
    Write(*,*)'Number of wavenumber bins = ', nkbins
    
    kindex_max = int((kmax-kmin)/dkappa)
    !kindex_max = int(kmax/dkappa)
    print *,'k0,kmax,kindex_max',k0,kmax, kindex_max!,eta, dchar/eta, kmax/(dchar/eta)
    
    !allocate(E_spectra(nkbins), karray(nkbins))
    
    allocate(E_spectra(0:kindex_max+1),karray(0:kindex_max+1))

    karray = zero 
    E_spectra  = zero 

    !print *,size(wx,1), size(wy,1), size(wz,1), size(E_spectra,1)
    
!!$    do k = 1, mz
!!$       do j = 1,my 
!!$          do i = 0, mx1
!!$             if(fluid_atijk(i,j,k)) then
!!$                particle_pos(i,j,k) = 1
!!$             else
!!$                particle_pos(i,j,k) = 0
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$    Write(*,'(A25,4(2x,g12.5))')'SUPERFICIAL VELOCITY = ', UMEAN
    
!!$    umean_int = zero 
!!$    count = 0
!!$    do i=1, mx1
!!$       do k=1,mz
!!$          do j=1,my
!!$             ur_tot(i,j,k,:) = ur_tot(i,j,k,:) + umean(:)
!!$             if (fluid_atijk(i,j,k))then 
!!$                count = count + 1
!!$                do idim=1,ndim
!!$                   umean_int(idim) = umean_int(idim) + ur_tot(i,j,k,idim)
!!$                end do
!!$             end if
!!$
!!$          end do
!!$       end do
!!$    end do

    !REMEMBER UMEAN IS THE SUPERFICIAL VELOCITY (OBTAINED BY AVERAGING OVER THE ENTIRE BOX VOLUME) AND UMEAN_INT IS THE INTERSTITIAL MEAN VELOCITY
!!$    umean_int(:) = umean_int(:)/count
!!$    Write(*,'(A25,3(2x,g12.5))') 'USF/(1-eps) = ', UMEAN/(one-vol_frac1)
!!$    Write(*,'(A25,3(2x,g12.5))')'UMEAN OF FLUCTUATING FIELD = ', UMEAN_INT
!!$    
!!$    Write(*,*)'Meanslipmod = ', meanslipmod

!!$    count = 0 
!!$    tke_r = zero 

!!$       do k=1,mz
!!$          do j=1,my
!!$             do i=1,mx1
!!$                ur_tot(i,j,k,idim) = (ur_tot(i,j,k,idim)) &
!!$                     & * particle_pos(i,j,k)
!!$                !ur_tot(i,j,k,idim) = (ur_tot(i,j,k,idim)-umean_int(idim)) &
!!$                !    & * particle_pos(i,j,k)
!!$                
!!$                arr(i,j,k,idim) = ur_tot(i,j,k,idim)!/ucharmod
!!$                
!!$                tke_r(idim) = tke_r(idim) + arr(i,j,k,idim)**2.d0
!!$
!!$             end do
!!$          end do
!!$       end do
    
    !turb_uf(:,:,:,1:3) = turb_uf(:,:,:,1:3)*(real(mx1*my*mz,prcn)/real(count_fluid,prcn))
    !count_real = count_fluid!SUM(particle_pos(1:mx1,1:my,1:mz))
!!$    count_real = mx1*my*mz
!!$    write(*,'(A25,5(2x,g12.7))') ' tke real',tke_r/count_real, sum(tke_r(:))/count_real
    
!!$    tke = 0.d0
!!$    do k=1,mz
!!$       do j=1,my
!!$          do i=1,my/2+1
!!$             if (i==1) then
!!$                tke(1) = tke(1)+dble(turb_uf(i,j,k,1)*conjg(turb_uf(i,j,k,1)))
!!$                tke(2) = tke(2)+dble(turb_uf(i,j,k,2)*conjg(turb_uf(i,j,k,2)))
!!$                tke(3) = tke(3)+dble(turb_uf(i,j,k,3)*conjg(turb_uf(i,j,k,3)))
!!$             else
!!$                tke(1) = tke(1)+dble(turb_uf(i,j,k,1)*conjg(turb_uf(i,j,k,1))) *2.d0
!!$                tke(2) = tke(2)+dble(turb_uf(i,j,k,2)*conjg(turb_uf(i,j,k,2))) *2.d0
!!$                tke(3) = tke(3)+dble(turb_uf(i,j,k,3)*conjg(turb_uf(i,j,k,3))) *2.d0
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$    count_comp =1! count_fluid!(mx1*my*mz)**2.d0!1!count_fluid
!!$    write(*,'(A25,5(2x,g12.7))') ' tke complex',tke/dble(count_comp), sum(tke(:))/dble(count_comp)
!!$    !write(*,'(A,4D20.12)') 'tke inside fixed bed(complex)',tke(:),sum(tke(:))
!!$
!!$    write(*,*) 'difference',sum(tke(:))/(sum(tke_r(:))/(count_real))
!!$    !  write(*,*)  size(E_spectra,1)

    do idim=1,ndim
       call ff3rc(ur_tot(1:mx1,1:my,1:mz,idim), turb_uf(1:my2,1:my,1:mz,idim))
    end do

    E_spectra = 0.d0

    epsf = 0.d0
    tketotal = 0.d0
    kmode_max = -1.0
    
    
    do k=1,mz
       do j=1,my
          do i=1,my/2+1
             tmp = -dble(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k))
             kmode = sqrt(tmp)

             if(kmode.eq.zero) goto 101 ! Corresponds to the mean. No
             ! flucutating energy in the mean!!
             
             if(kmode.gt.kmode_max) kmode_max = kmode
             if (kmode >= kmax) then
                kmode = kmax
             end if
             kindex = int((kmode-kmin)/dkappa)

             karray(int(kindex)) = kmin + kindex*dkappa

             tmp = 0.d0

             do idim=1,ndim
                tmp = tmp+dble(turb_uf(i,j,k,idim)*conjg(turb_uf(i,j,k,idim)))
             end do
             
             if (i==1) then
                E_spectra(int(kindex)) = E_spectra(int(kindex)) + tmp
                
                tketotal = tketotal + tmp
             else
                E_spectra(int(kindex)) = E_spectra(int(kindex)) + tmp*2.d0
                
                tketotal = tketotal + tmp*2.d0
             end if
101          CONTINUE
          end do
       end do
    end do
    !  print *,'after calculate E_spectra'

    Write(*,*) ' Maximum KMODE = ', kmode_max
    
    write(*,'(A25,5(2x,g12.7))') ' tke from spectra ',sum(E_spectra(:))/dble(count_comp)
    !print *,'epsf=',epsf
    open(30,file=TRIM(RUN_NAME)//'_E_spectra.dat',form='formatted')
    
    write(30 , *)  'zone t=" total bed  " '
    
    
!!$    lenscale = ((1-maxvolfrac)*lx*ly*lz)**(1.d0/3.d0)
!!$    !lenscale = (nbody*pi*dchar*dchar*dchar/6.d0)**(1.d0/3.d0) 
!!$    write(*,'(A25,g17.8)') 'Number of bodies = ', real(nbody)
!!$    write(*,'(A25,g17.8)') 'Normalizing by total tke = ', tketotal
!!$    write(*,'(A25,g17.8)') 'length scale defined as (\eps * V)^1/3 = ', lenscale
!!$    write(30,'(A40)') '#kappa,kappa*lx,kappa*d,kappa*len,E,E/tketotal'
!!$    write(*,'(A25,g17.8)') 'Diameter of the particle = ', dchar
!!$    write(*,'(A60,g17.8, I)') 'Wavenumber corresponding to the diamete&
!!$         &r of the particle = ', two*pi/dchar, int(two*pi/(dchar&
!!$         &*dkappa))
    
    Write(*,*) ' imis = ', imis
    
    do i=0, kindex_max
       
       E_spectrum_function = E_spectra(i)
       
       write(30,'(10(D20.12))') karray(i), karray(i)*dchar/(two*pi),&
            & E_spectrum_function, karray(i)**(-5.d0&
            &/3.d0)
       
       karray_mis(i,imis) = karray(i)
       E_spectrum_function_mis(i,imis) = E_spectrum_function/(dkappa&
            &*tke_mis(imis)*dchar/(two*pi))!mean_energy
    end do
    nkbins_mis = nkbins

    close(30,status="keep")
    Write(*,*) 'I AM HERE NOW'
    !deallocate(E_spectra,karray, kcount)
    
    !if(iscalon.eq.1) CALL calc_scalar_spec_3d
  end SUBROUTINE calc_spec_3d
  
  
  SUBROUTINE calc_scalar_spec_3d
    implicit none
    double complex, dimension(my2,my,mz,nspmx) :: turb_phif
    double complex, dimension(my2,mz) :: phif_tmp
    double precision, dimension(my,mz) :: phir_tmp
    double complex :: wx(my2), wy(my), wz(mz)
    double precision, allocatable, dimension(:) :: karray
    double precision, allocatable, dimension(:,:) :: E_spectra
    !double precision, allocatable, dimension(:,:,:) :: arr
    integer :: particle_pos(0:mx+1,my,mz)
    double precision :: k0, kmax, tmp(nspmx), kmode, lx, ly, lz, tmp1, lenscale, tketotal(nspmx)
    double precision :: phimean_int(nspmx), phi_tke(nspmx),  phi_tker(nspmx)
    double precision :: kmaxeta, eta, dkappa
    integer :: i,j,k,isp, count,turb_lx, ibody, istat, count_real, count_comp,kindex, kindex_max
    
    integer :: iloc, jloc, kloc, xlo, xhi, ylo, yhi, zlo, zhi, inpos, outpos
    double precision :: dist, energy_bal, epsf, ktot


    ly = doml(2)
    lz = doml(3)
    lx = doml(1)
    Write(*,*) 'Computing scalar spectra'
    print *,'lx,ly,lz',lx,ly,lz

    do i=1,my
       do isp=1,nspmx
          call ff2cr(phif(i,:,:,isp), ur_tot(i,:,:,isp))
       end do
    end do

    
    do i=1, my/2+1
       wx(i)=dcmplx(zero,2.0*pi*dble(i-1)/dble(lx))
    enddo
    
    do i=1,my/2
       wy(i) = dcmplx(zero,2*pi*dble(i-1)/dble(ly))
       wy(my+1-i) = dcmplx(zero,-1.d0*2*pi*dble(i)/dble(ly))
    enddo
    
    do i=1,mz/2
       wz(i)=dcmplx(zero,2.0*pi*dble(i-1)/dble(lz))
       wz(mz+1-i)=dcmplx(zero,-1.d0*2*pi*dble(i)/dble(lz))
    enddo
    
    
    kmaxeta = 1.0
    k0 = aimag(wx(2))
    !kmax = sqrt(2.0)/3.d0*my*k0
    kmax = my*k0
    eta = kmaxeta/kmax
    !dkappa shud be in the radial/spherical coordinate system. THe above wavenumbers are defined
    !on a cartesian grid and dkappa is chosen in the radial system
    dkappa = dsqrt(3.d0)*two*pi/lx
    kindex_max = int(kmax/dkappa)
    print *,'k0,kmax,eta,d/eta, k_d',k0,kmax, eta, dchar/eta, kmax/(dchar/eta)
    
    allocate(E_spectra(0:kindex_max+1,1:nspmx), karray(0:kindex_max+1))
    karray = zero 
    E_spectra  = zero 
    !print *,size(wx,1), size(wy,1), size(wz,1), size(E_spectra,1)
    
    do k = 1, mz
       do j = 1,my 
          do i = 0, mx1
             if(fluid_atijk(i,j,k)) then
                particle_pos(i,j,k) = 1
             else
                particle_pos(i,j,k) = 0
             end if
          end do
       end do
    end do
    
    
    phimean_int = zero 
    count = 0
    do i=1, mx1
       do k=1,mz
          do j=1,my
             ur_tot(i,j,k,1:nspmx) = ur_tot(i,j,k,1:nspmx) + phirmean(1:nspmx)
             if (fluid_atijk(i,j,k))then 
                count = count + 1
                do isp=1,nspmx
                   phimean_int(isp) = phimean_int(isp) + ur_tot(i,j,k,isp)
                end do
             end if
             
          end do
       end do
    end do

    phimean_int(:) = phimean_int(:)/real(count,prcn)
    
    
    count = 0 

    do isp=1,nspmx
       do k=1,mz
          do j=1,my
             do i=1,mx1
                ur_tot(i,j,k,isp) = (ur_tot(i,j,k,isp)-phimean_int(isp)) &
                     & * particle_pos(i,j,k)
                
                arr(i,j,k,isp) = ur_tot(i,j,k,isp)
                               
             end do
          end do
       end do
       call ff3rc(arr(1:mx1,1:my,1:mz,isp), turb_phif(1:my2,1:my,1:mz,isp))
       
    end do

    E_spectra = 0.d0
    
    epsf = 0.d0
    tketotal = 0.d0
    do k=1,mz
       do j=1,my
          do i=1,my/2+1
             tmp(1) = -dble(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k))
             kmode = sqrt(tmp(1))
             if (kmode >= kmax) then
                kmode = kmax
             end if
             kindex = int(kmode/dkappa)
             karray(int(kindex)) = kindex*dkappa
             !PRINT*,'kmode = ', kmode, sqrt(tmp)
             tmp = 0.d0
             do isp=1,nspmx
                tmp(isp) = dble(turb_phif(i,j,k,isp)*conjg(turb_phif(i,j,k,isp)))
                
                if (i==1) then
                   E_spectra(int(kindex),isp) = E_spectra(int(kindex),isp) + tmp(isp)
                   
                   tketotal(isp) = tketotal(isp) + tmp(isp)
                else
                   E_spectra(int(kindex),isp) = E_spectra(int(kindex),isp) + tmp(isp)*2.d0
                   tketotal(isp) = tketotal(isp) + tmp(isp)*2.d0
                end if
             end do
          end do
       end do
    end do
    !  print *,'after calculate E_spectra'

    open(30,file=TRIM(POST_RUNNAME)//'_PHI_spectra.dat',form='formatted')

    write(30 , *)  'zone t=" total bed  " '
    
    ktot = zero
    lenscale = ((1-maxvolfrac)*lx*ly*lz)**(1.d0/3.d0)
    !lenscale = (nbody*pi*dchar*dchar*dchar/6.d0)**(1.d0/3.d0) 
    write(*,'(A25,g17.8)') 'Number of bodies = ', real(nbody)
    write(*,'(A25,g17.8)') 'Normalizing by total tke = ', tketotal
    write(*,'(A25,g17.8)') 'length scale defined as (\eps * V)^1/3 = ', lenscale
    write(30,'(A40)') '#kappa,kappa*lx,kappa*d,kappa*len,E,E/tketotal'
    do i=0, kindex_max
       !PRINT*,'karrayt = ', karray(i)
       if(E_spectra(i,1).gt.zero) then
          write(30,'(10(D20.12))') karray(i), karray(i)*dchar/(maxvolfrac), karray(i)*dchar/(two*pi), karray(i)*lenscale, E_spectra(i,1), E_spectra(i,1)/tketotal(1) 
       end if
    end do
    ! WRITE(*,*)'ktot = ', ktot
    close(30,status="keep")
    deallocate(E_spectra, karray)
  end SUBROUTINE calc_scalar_spec_3d


  SUBROUTINE calc_fluid_scalar_stats
    USE dependent_functions
    USE nlarrays
    USE lang_source_diss
    IMPLICIT NONE

    REAL(prcn) :: umean_int(ndim),ustdev(ndim), tke, sdevu,&
         & scalar_var(nspmx), phimean_int(nspmx), velvar(ndim)
    INTEGER :: i, j, k, idim, isp, someunit, flag_uphi(ndim), jdim
    REAL(prcn) :: scalar_vel_corr(nspmx,ndim), phistdev(nspmx), umag,&
         & umean_intmod, fluidfluctvel(ndim), tke_par, tke_per,&
         & uparallel, uperp(ndim), uiuj(ndim,ndim), xsi, eta

    CALL calc_velreal(ur_tot)

    if(iscalon.eq.1)then
       do isp=1,nspmx
          do i=1,mx1,1
             do j=1,my2,1
                do k=1,mz,1
                   uf1(j,k) = phif(i,j,k,isp)
                enddo
             enddo

             call ff2cr(uf1,ur1)	 
             do j=1,my,1
                do k=1,mz,1
                   arr(i,j,k,isp) = ur1(j,k)+phirmean(isp)
                enddo
             enddo
          enddo
       enddo
    end if

    umean_int = zero
    ustdev = zero

    phimean_int = zero
    phistdev = zero

    do i=1, mx1
       do k=1,mz
          do j=1,my
             if (fluid_atijk(i,j,k))then 
                do idim=1,ndim
                   umean_int(idim) = umean_int(idim) + ur_tot(i,j,k&
                        &,idim)
                   ustdev(idim) = ustdev(idim) + ur_tot(i,j,k,idim)&
                        &**two
                end do
                if(iscalon.eq.1)then
                   do isp = 1, nspmx
                      phimean_int(isp) = phimean_int(isp) + arr(i,j,k&
                           &,isp)
                      phistdev(isp) = phistdev(isp) + arr(i,j,k,isp)&
                           &**two
                   end do
                end if
             end if
          end do
       end do
    end do
    umean_int(:) = umean_int(:)/real(count_fluid,prcn)
    ustdev(:) = DSQRT(ustdev(:)/real(count_fluid,prcn) - umean_int(:)&
         &**two)
    
    umean_intmod = DOT_PRODUCT(umean_int(1:ndim),umean_int(1:ndim))
    umean_intmod = DSQRT(umean_intmod)
    
    parallel_dir(1:ndim) = umean_int(1:ndim)/umean_intmod
    
    if(iscalon.eq.1)then
       phimean_int(:) = phimean_int(:)/real(count_fluid,prcn)
       phistdev(:) = DSQRT(phistdev(:)/real(count_fluid,prcn) -&
            & phimean_int(:)**two)
    end if
    !    sdevu = DOT_PRODUCT(ustdev(1:ndim),ustdev(1:ndim))
    sdevu = DOT_PRODUCT(umean_int(1:ndim),umean_int(1:ndim))
    sdevu = DSQRT(sdevu)
    tke = zero 
    tke_par = zero
    tke_per = zero
    velvar = zero
    scalar_var = zero
    scalar_vel_corr = zero
    uiuj = zero
    
    do k=1,mz
       do j=1,my
          do i=1,mx1
             if(fluid_atijk(i,j,k))then 
                
                do idim = 1, ndim
                   velvar(idim) = velvar(idim) + (ur_tot(i,j,k,idim)&
                        &-umean_int(idim))**2.d0
                end do

                fluidfluctvel(1:ndim) = ur_tot(i,j,k,1:ndim)&
                     &-umean_int(1:ndim)

                uparallel = DOT_PRODUCT(fluidfluctvel(1:ndim)&
                     &,parallel_dir(1:ndim))
                
                tke_par = tke_par + uparallel**2.d0

                uperp(1:ndim) = fluidfluctvel(1:ndim) - uparallel&
                     &*parallel_dir(1:ndim)
                
                tke_per = tke_per + DOT_PRODUCT(uperp(1:ndim),&
                     & uperp(1:ndim))

                do idim = 1, ndim
                   do jdim = 1, ndim
                      uiuj(idim,jdim) = uiuj(idim,jdim) +&
                           & fluidfluctvel(idim)*fluidfluctvel(jdim)
                   end do
                end do
                
                ur_tot(i,j,k,1:ndim) = ur_tot(i,j,k,1:ndim) -&
                     & umean_int(1:ndim)
                
                if(iscalon.eq.1)then
                   do isp = 1, nspmx
                      scalar_var(isp) = scalar_var(isp)+(arr(i,j,k,isp)-phimean_int(isp))**2.d0
                      do idim = 1, ndim
                         scalar_vel_corr(isp,idim) = scalar_vel_corr(isp,idim)+(arr(i,j,k,isp)-phimean_int(isp))*(ur_tot(i,j,k,idim)-umean_int(idim))
                      end do
                   end do
                   
                   arr(i,j,k,1:nspmx) = arr(i,j,k,1:nspmx) - phimean_int(1:nspmx)
                end if
             else
                ur_tot(i,j,k,1:ndim) = zero
                if(iscalon.eq.1) arr(i,j,k,1:nspmx) = zero
             endif
          end do
       end do
    end do

    velvar(1:ndim) = (velvar(1:ndim)/real(count_fluid,prcn))
    
    do idim = 1, ndim
       do jdim = 1, ndim
          uiuj(idim,jdim) = uiuj(idim,jdim)/real(count_fluid,prcn)
       end do
    end do

    if(iscalon.eq.1) then
       scalar_var(1:nspmx) = scalar_var(1:nspmx)/real(count_fluid,prcn)
       do isp = 1, nspmx
          do idim = 1, ndim
             scalar_vel_corr(isp, idim) = scalar_vel_corr(isp, idim)/real(count_fluid,prcn)
             scalar_vel_corr(isp, idim) = scalar_vel_corr(isp, idim)/DSQRT(velvar(idim)*scalar_var(isp))
          end do
       end do
    end if

    tke = (velvar(1) + velvar(2) + velvar(3))

    tke_par = tke_par/real(count_fluid,prcn)
    tke_per = tke_per/real(count_fluid,prcn)

    someunit = getnewunit(minunitno,maxunitno)
    
    OPEN(someunit, FILE=TRIM(RUN_NAME)//'_scal_stats.dat', form='formatted')
    Write(someunit, '(9(2x,g17.5))') velvar(1)/meanslipmod**2.d0, velvar(2)/meanslipmod**2.d0, velvar(3)/meanslipmod**2.d0, tke/meanslipmod**2.d0, scalar_var(1), (scalar_vel_corr(1,idim), idim = 1, ndim)
    CLOSE(someunit, status='keep')
    
    tke_mis(imis) = tke/umean_intmod**2.d0!ucharmod**2.d0

    tke_par_mis(imis) = tke_par/umean_intmod**2.d0

    tke_per_mis(imis) = tke_per/umean_intmod**2.d0

    scalar_var_mis(imis) = scalar_var(1)
    scalar_vel_corr_mis(1:ndim,imis) = scalar_vel_corr(1,1:ndim)


    CALL calc_anisotropy(uiuj(1:ndim,1:ndim), xsi, eta)

    xsi_mis(imis) = xsi
    eta_mis(imis) = eta

!!$    OPEN(someunit, FILE=TRIM(POST_RUNNAME)//'_u_phi_fluct_fields.dat', form='formatted')
!!$    write(someunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "phi_n" ' ,' &
!!$         &"UX" ',' "UY" ',' "UZ" ' ,' "FLAGX" ' ,' "FLAGY" ' ,' "FLAGZ" ' 
!!$    write(someunit,*)'ZONE F=POINT, I=', mx1,  ', J=', my, ', K=', mz
!!$    do k = 1, mz
!!$       do j = 1, my
!!$          do i = 1, mx1
!!$             flag_uphi(1:ndim) = 0
!!$             do idim = 1, ndim
!!$                if((ur_tot(i,j,k,idim)*arr(i,j,k,1)).gt.zero) then
!!$                   flag_uphi(idim) = 1
!!$                else if((ur_tot(i,j,k,idim)*arr(i,j,k,1)).lt.zero) then
!!$                   flag_uphi(idim) = -1
!!$                end if
!!$              end do
!!$              umag = DSQRT(DOT_PRODUCT(ur_tot(i,j,k,1:ndim), ur_tot(i,j,k,1:ndim)))
!!$                write(someunit,*) real((i-1))*dx/dchar,real((j-1))*dx/dchar, (k&
!!$                     &-1)*dx/dchar,  arr(i,j,k,1)/DSQRT(scalar_var(1)), ur_tot(i,j,k,1)/DSQRT(velvar(1)) ,ur_tot(i,j,k,2)/DSQRT(velvar(2))&
!!$                     &,ur_tot(i,j,k,3)/DSQRT(velvar(3)) , umag/DSQRT(tke), flag_uphi(1), flag_uphi(2), flag_uphi(3)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$21  FORMAT(10(1xe17.4))
!!$    CLOSE(someunit, status='keep')

    RETURN

    Write(*,*) ' Computing the Energy Spectrum '

    CALL calc_spec_3d

    Write(*,*) ' Done with computing the Energy Spectrum '

    Write(*,*) 'Computing fluid velocity PDFs'

    CALL calc_fluid_vel_pdfs

    Write(*,*) ' Done with computing fluid velocity PDFs'
  END SUBROUTINE calc_fluid_scalar_stats

  SUBROUTINE calc_mis_avg_vel_field
    USE dependent_functions

    IMPLICIT NONE
    
    INTEGER :: i, j, k
    
    CALL calc_velreal(ur_tot)
    
    do k = 1, mz
       do j = 1, my
          do i = 1, mx1
             if(fluid_atijk(i,j,k))then
                misvel(i,j,k,1:ndim) = misvel(i,j,k,1:ndim) + ur_tot(i,j,k,1:ndim)
             else
                misvel(i,j,k,1:ndim) = misvel(i,j,k,1:ndim) + zero
             end if
          end do
       end do
    end do
    
  END SUBROUTINE calc_mis_avg_vel_field

  SUBROUTINE calc_filtered_mis_avg_fields
    IMPLICIT NONE
    
    INTEGER :: i, j, k, idim, misavgenerunit
    REAL(prcn) ::  varvel(ndim), tke
    
    misavg_mean_vel(1:ndim) = zero

    do idim = 1, ndim
       do k = 1, mz
          do j = 1, my
             do i = 1, mx1
                misvel(i,j,k,idim) = misvel(i,j,k,idim)/(real(countmis,prcn)*(one-vol_frac1))
                misavg_mean_vel(idim) = misavg_mean_vel(idim) + misvel(i,j,k,idim)
             end do
          end do
       end do
    end do

    misavg_mean_vel(1:ndim) = misavg_mean_vel(1:ndim)/real(mx1*my*mz, prcn)
    
    Write(*,*)' MIS AVG MEAN VELOCITY : ', misavg_mean_vel(1), ucharmod/(one-vol_frac1)
    
    CALL calc_mis_avg_field_spectra
!!$    varvel(1:ndim) = zero
!!$    do idim = 1, ndim
!!$       do k = 1, mz
!!$          do j = 1, my
!!$             do i = 1, mx1
!!$                varvel(idim) = varvel(idim) + (arr(i,j,k,idim)-misavg_mean_vel(idim))**2.d0
!!$             end do
!!$          end do
!!$       end do
!!$       varvel(idim) = varvel(idim)/real(mx1*my*mz, prcn)
!!$    end do
!!$
!!$    tke = SUM(varvel(1:ndim))
!!$    
!!$    misavgenerunit = getnewunit(minunitno,maxunitno)
!!$
!!$    OPEN(misavgenerunit, FILE=TRIM(POST_RUNNAME)//'_MISAVG_ENER.dat', form='formatted')
!!$    Write(misavgenerunit, '(9(2x,g17.5))') varvel(1)/ucharmod**2.d0, varvel(2)/ucharmod**2.d0, varvel(3)/ucharmod**2.d0, tke/ucharmod**2.d0
!!$    CLOSE(misavgenerunit, status='keep')
    
  END SUBROUTINE calc_filtered_mis_avg_fields

  SUBROUTINE calc_mis_avg_field_spectra
    implicit none
    double complex, dimension(my2,my,mz,ndim) :: turb_uf
    double complex :: wx(my2), wy(my), wz(mz)
    double precision, allocatable, dimension(:) :: E_spectra, karray
    double precision :: k0, kmax, tmp, kmode, lx, ly, lz, tmp1, lenscale, tketotal, kmin 
    double precision :: tke(ndim),  tke_r(ndim)
    double precision :: kmaxeta, eta, dkappa
    integer :: i,j,k,idim, count,turb_lx, ibody, istat, count_real, count_comp,kindex, kindex_max

    integer :: iloc, jloc, kloc, xlo, xhi, ylo, yhi, zlo, zhi, inpos, outpos, nkbins
    double precision :: dist, energy_bal, epsf, ktot, kmode_max,&
         & dxlength, avg_phii, decimal, E_spectrum_function

    REAL(prcn) :: mean_energy

    ly = doml(2)
    lz = doml(3)
    lx = doml(1)
    print *,'lx,ly,lz',lx,ly,lz

    mean_energy = (ucharmod/(one-vol_frac1))**2.d0

    do i=1, my/2+1
       wx(i)=dcmplx(zero,2.0*pi*dble(i-1)/dble(lx))
    enddo
    
    do i=1,my/2
       wy(i) = dcmplx(zero,2*pi*dble(i-1)/dble(ly))
       wy(my+1-i) = dcmplx(zero,-1.d0*2*pi*dble(i)/dble(ly))
    enddo
    
    do i=1,mz/2
       wz(i)=dcmplx(zero,2.0*pi*dble(i-1)/dble(lz))
       wz(mz+1-i)=dcmplx(zero,-1.d0*2*pi*dble(i)/dble(lz))
    enddo
    
    k0 = aimag(wx(2))
    !kmax = sqrt(2.0)/3.d0*my*k0
    kmax = my*k0/two ! Smallest wavelength that can be resolved is of grid size 2*dx (not dx). 
    !dkappa shud be in the radial/spherical coordinate system. THe
    ! above wavenumbers are defined on a cartesian grid and
    ! dkappa is chosen in the radial system

    dkappa = dsqrt(3.d0)*two*pi/lx
    Write(*,*)'dkappa = ', dkappa
    
    
    !nkbins = INT((kmax-kmin)/dkappa)
    nkbins = INT(kmax/dkappa)
    Write(*,*)'Number of wavenumber bins = ', nkbins
    
    
    kindex_max = int(kmax/dkappa)
    print *,'k0,kmax,kindex_max',k0,kmax, kindex_max!,eta, dchar/eta, kmax/(dchar/eta)
    allocate(E_spectra(0:kindex_max+1), karray(0:kindex_max+1))
    
    karray = zero 
    E_spectra  = zero 

    count = 0 
    tke_r = zero 
    do idim=1,ndim
       do k=1,mz
          do j=1,my
             do i=1,mx1
!!$                misvel(i,j,k,idim) = misvel(i,j,k,idim)&
!!$                     &-misavg_mean_vel(idim)

                misvel(i,j,k,idim) = misvel(i,j,k,idim)

                misvel(i,j,k,idim) = misvel(i,j,k,idim)!/ucharmod*&
                !& (one-vol_frac1)
             end do
          end do
       end do
       call ff3rc(misvel(1:mx1,1:my,1:mz,idim), turb_uf(1:my2,1:my,1:mz,idim))

    end do

    E_spectra = 0.d0
    
    do k=1,mz
       do j=1,my
          do i=1,my/2+1
             tmp = -dble(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k))
             kmode = sqrt(tmp)
             !if(kmode.eq.zero) goto 101 ! Corresponds to the mean. No
             ! flucutating energy in the mean!!
             
             if(kmode.gt.kmode_max) kmode_max = kmode
             if (kmode >= kmax) then
                kmode = kmax
             end if
             kindex = int(kmode/dkappa)
             karray(int(kindex)) = kindex*dkappa
             
             tmp = 0.d0
             do idim=1,ndim
                tmp = tmp+dble(turb_uf(i,j,k,idim)*conjg(turb_uf(i,j,k,idim)))
             end do
             
             if (i==1) then
                E_spectra(int(kindex)) = E_spectra(int(kindex)) + tmp
                
                tketotal = tketotal + tmp
             else
                E_spectra(int(kindex)) = E_spectra(int(kindex)) + tmp*2.d0
                
                tketotal = tketotal + tmp*2.d0
             end if
101          CONTINUE
          end do
       end do
    end do

    
    Write(*,*) ' Maximum KMODE = ', kmode_max
    
    !print *,'epsf=',epsf
    open(30,file=TRIM(POST_RUNNAME)//'_MISAVG_FIELD_spectra.dat',form='formatted')
    
    write(30 , *)  'zone t=" total bed  " '
    
    do i=0, kindex_max
       
       E_spectrum_function = E_spectra(i)
       
       write(30,'(10(D20.12))') karray(i), karray(i)*dchar/(two*pi),&
            & E_spectrum_function/mean_energy, karray(i)**(-5.d0&
            &/3.d0)
    end do
    
    close(30,status="keep")
    
    deallocate(E_spectra,karray)
    
  end SUBROUTINE calc_mis_avg_field_spectra


  SUBROUTINE calc_residual_energy
    IMPLICIT NONE
    
    INTEGER :: imeasvol, idim, measvolcount, ICELLS, JCELLS, KCELLS, i, j, k
    REAL(prcn) :: Lmeas, Lmeasmin, Lmeasmax

    CALL calc_velreal(ur_tot)

    do k = 1, mz
       do j = 1, my
          do i = 1, mx1
             if(fluid_atijk(i,j,k))then
                meanuf(1:ndim) = meanuf(1:ndim) + ur_tot(i,j,k,1:ndim)
             end if
          end do
       end do
    end do
    meanuf(1:ndim) = meanuf(1:ndim)/real(count_fluid,prcn)

    Lmeasmax = DOML(2)
    Lmeasmin = dia_phys
    Lmeas = Lmeasmax
    
    measvolcount = 0
    
    Do While(Lmeas.gt.Lmeasmin)
       
       IF(ASSOCIATED(COARSE_GRID%XE))DEALLOCATE(COARSE_GRID%XE)
       IF(ASSOCIATED(COARSE_GRID%YN))DEALLOCATE(COARSE_GRID%YN)
       IF(ASSOCIATED(COARSE_GRID%ZT))DEALLOCATE(COARSE_GRID%ZT)
       
       CALL INITIALIZE_GRID(Lmeas, ICELLS, JCELLS, KCELLS)
       
       CALL COMPUTE_ENERGY_IN_FILTER_WIDTH(Lmeas) 

       measvolcount = measvolcount + 1       

       energy_measvol(imis,measvolcount)= COARSE_GRID%ener 
       energy_res_measvol(imis,measvolcount) = COARSE_GRID%ener_res
       energy_cross_measvol(imis,measvolcount) = COARSE_GRID%ener_cross

       Lmeasvol(measvolcount) = Lmeas

       Lmeas = DOML(2)/real(measvolcount+1,prcn)
       
    End Do
    
    totalmeasvolcount = measvolcount

  END SUBROUTINE calc_residual_energy
  
  SUBROUTINE compute_energy_in_filter_width(Lmeas)
    IMPLICIT NONE
    REAL(prcn), Intent(in) :: Lmeas
    INTEGER :: cx, cy, cz, i, imin, imax, jmin, jmax, kmin, kmax, ifinew, ifinee, jfines, jfinen, kfineb, kfinet, j, k, ii, jj, kk, PIJK(ndim), m, idim, count_fluid_measvol, count_coarse
    REAL(prcn) :: xe, xw, ys, yn, zb, zt
    REAL(prcn), ALLOCATABLE, DIMENSION(:,:,:,:) :: meanuf_measvol, fluctuf_measvol
    REAL(prcn), ALLOCATABLE, DIMENSION(:,:,:) :: ener_measvol, ener_res_measvol, ener_cross_measvol
    LOGICAL :: ALREADY_COUNTED(mx1,my,mz)
    REAL(prcn) :: mean_enerf, ener, ener_res, ener_cross
    
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

    ALLOCATE(meanuf_measvol(IMAX,JMAX,KMAX,NDIM), fluctuf_measvol(IMAX, JMAX, KMAX,NDIM))
    ALLOCATE(ener_measvol(IMAX,JMAX,KMAX),ener_res_measvol(IMAX,JMAX,KMAX), ener_cross_measvol(IMAX,JMAX,KMAX))

    do k = KMIN, KMAX
       ZT = COARSE_GRID%ZT(K)
       ZB = ZT - COARSE_GRID%dz
       !Write(*,*)'ZT = ', ZT, ZB
       KFINEB = ceiling(ZB/dz) + 1
       KFINET = floor(ZT/dz) + 1
       !Write(*,*) ' KFINES = ', KFINEB, KFINET
       do j = JMIN, JMAX
          YN = COARSE_GRID%YN(J)
          YS = YN - COARSE_GRID%dy
          JFINES = ceiling(YS/dy) + 1
          JFINEN = floor(YN/dy) + 1
          !Write(*,*) ' JFINES = ', JFINES, JFINEN
          do i = IMIN, IMAX
             XE = COARSE_GRID%XE(I)
             XW = XE - COARSE_GRID%dx
             IFINEW = ceiling(XW/dx) + 1
             IFINEE = floor(XE/dx) + 1
             !Write(*,*) ' IFINES =', IFINEW, IFINEE
             count_fluid_measvol = 0 
             meanuf_measvol(i,j,k,1:ndim) = zero
             do KK = KFINEB, KFINET
                do JJ = JFINES, JFINEN
                   do II = IFINEW, IFINEE
                      if((KK.le.mz).and.(JJ.le.my).and.(II.le.mx1))then
                         if(fluid_atijk(II,JJ,KK))then
                            if(.not.ALREADY_COUNTED(II,JJ,KK))then
                               meanuf_measvol(i,j,k,1:ndim) = meanuf_measvol(i,j,k,1:ndim) + ur_tot(II,JJ,KK,1:ndim)
                               count_fluid_measvol = count_fluid_measvol + 1
                               ALREADY_COUNTED(II,JJ,KK) = .TRUE.
                            end if
                         end if
                      end if
                   end do
                end do
             end do
             meanuf_measvol(i,j,k,1:ndim) = meanuf_measvol(i,j,k,1:ndim)/real(count_fluid_measvol,prcn)
             fluctuf_measvol(i,j,k,1:ndim) = meanuf(1:ndim) - meanuf_measvol(i,j,k,1:ndim)
          end do
       end do
    end do
    mean_enerf = half*DOT_PRODUCT(meanuf(1:ndim), meanuf(1:ndim))

    count_coarse = 0
    ener = zero
    ener_res = zero
    ener_cross = zero
    do k = kmin, kmax
       do j = jmin, jmax
          do i = imin, imax
             ener_measvol(i,j,k) = half*DOT_PRODUCT(meanuf_measvol(i,j,k,1:ndim), meanuf_measvol(i,j,k,1:ndim))
             ener_res_measvol(i,j,k) = half*DOT_PRODUCT(fluctuf_measvol(i,j,k,1:ndim), fluctuf_measvol(i,j,k,1:ndim))
             ener_cross_measvol(i,j,k) = DOT_PRODUCT(meanuf_measvol(i,j,k,1:ndim), fluctuf_measvol(i,j,k,1:ndim))

             ener = ener + ener_measvol(i,j,k)
             ener_res = ener_res + ener_res_measvol(i,j,k)
             ener_cross = ener_cross + ener_cross_measvol(i,j,k)
             
             count_coarse = count_coarse + 1

!!$             Write(*,*)'Energy in ', i, j, k, 'is :', ener_measvol(i&
!!$                  &,j,k) + ener_res_measvol(i,j,k) +&
!!$                  & ener_cross_measvol(i,j,k), 'Total energy = ',&
!!$                  & mean_enerf
!!$
!!$             Write(*,*)'Individual Components of Energy in ', i, j, k&
!!$                  &, 'are :', ener_measvol(i,j,k), &
!!$                  & ener_res_measvol(i,j,k), ener_cross_measvol(i,j,k)
             
          end do
       end do
    end do
    ener = ener/real(count_coarse,prcn)
    ener_res = ener_res/real(count_coarse,prcn)
    ener_cross = ener_cross/real(count_coarse,prcn)

    COARSE_GRID%ener = ener/mean_enerf
    COARSE_GRID%ener_res = ener_res/mean_enerf
    COARSE_GRID%ener_cross = ener_cross/mean_enerf


    Write(*,*)'Energy in MEAS Length :', lmeas, 'is :',ener, ener_res , ener_cross
    
  END SUBROUTINE compute_energy_in_filter_width


  SUBROUTINE initialize_grid(Lmeas, ICELLS, JCELLS, KCELLS)
    IMPLICIT NONE
    REAL(prcn), Intent(in) :: Lmeas
    INTEGER, INTENT(out) :: ICELLS, JCELLS, KCELLS
    INTEGER :: cx, cy, cz, i, imin, jmin, kmin, ifinew,ifinee, jfines, jfinen, kfineb, kfinet, j, k

    cx = MAX(INT(DOML(1)/Lmeas)+1,2)
    cy = MAX(INT(DOML(2)/Lmeas)+1,2)
    cz = MAX(INT(DOML(3)/Lmeas)+1,2)

    COARSE_GRID%dx = DOML(1)/REAL(cx-1, prcn)
    COARSE_GRID%dy = DOML(2)/REAL(cy-1, prcn)
    COARSE_GRID%dz = DOML(3)/REAL(cz-1, prcn)

    COARSE_GRID%cx = cx
    COARSE_GRID%cy = cy
    COARSE_GRID%cz = cz

    ALLOCATE(coarse_grid%XE(cx-1),coarse_grid%YN(cy-1),coarse_grid%ZT(cz-1))

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

end MODULE spectra_3d
