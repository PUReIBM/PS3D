MODULE spectra_3d
  USE precision  
  USE constants  
  USE global_data
  USE nlmainarrays, ur_tot=>ubc, arr=>nlbc
  use fftw_interface

  implicit none
CONTAINS 
  SUBROUTINE calc_spec_3d
    double complex, dimension(my2,my,mz,ndim) :: turb_uf
    double complex, dimension(my2,mz) :: uf_tmp
    double precision, dimension(my,mz) :: ur_tmp
    double complex :: wx(my2), wy(my), wz(mz)
    !double precision, dimension(mx,my,mz,ndim) :: ur_tot
    double precision, allocatable, dimension(:) :: E_spectra, E_uu, karray
    !double precision, allocatable, dimension(:,:,:) :: arr
    integer :: particle_pos(0:mx+1,my,mz)
    double precision :: k0, kmax, tmp, kmode, lx, ly, lz, tmp1
    double precision :: umean_int(ndim), tke(ndim),  tke_r(ndim)
    double precision :: kmaxeta, eta, dkappa
    integer :: i,j,k,idim, count,foffset, nbody,turb_lx, ibody, istat, count_real, count_comp,kindex, kindex_max
  
    integer :: iloc, jloc, kloc, xlo, xhi, ylo, yhi, zlo, zhi, inpos, outpos
    double precision :: dist, energy_bal, epsf, ktot
    

    ly = doml(2)
    lz = doml(3)
    lx = doml(1)
    print *,'lx,ly,lz',lx,ly,lz

    do i=1,my
       do idim=1,ndim
          call ff2cr(u(i,:,:,idim), ur_tot(i,:,:,idim))
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
    kmax = sqrt(2.0)/3.d0*my*k0
    eta = kmaxeta/kmax
    dkappa = two*pi/lx
    kindex_max = int(kmax/dkappa)
    print *,'k0,kmax,eta,d/eta, k_d',k0,kmax, eta, dchar/eta, kmax/(dchar/eta)
    
    allocate(E_spectra(0:kindex_max+1), E_uu(0:kindex_max+1),karray(0:kindex_max+1))
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
    Write(*,'(A25,4(2x,g12.5))')'SUPERFICIAL VELOCITY = ', UMEAN

    umean_int = zero 
    count = 0
    do i=1, mx1
       do k=1,mz
          do j=1,my
             ur_tot(i,j,k,:) = ur_tot(i,j,k,:) + umean(:)
             if (fluid_atijk(i,j,k))then 
                count = count + 1
                do idim=1,ndim
                   umean_int(idim) = umean_int(idim) + ur_tot(i,j,k,idim)
                end do
             end if
             
          end do
       end do
    end do
    
    !REMEMBER UMEAN IS THE SUPERFICIAL VELOCITY (OBTAINED BY AVERAGING OVER THE ENTIRE BOX VOLUME) AND UMEAN_INT IS THE INTERSTITIAL MEAN VELOCITY
    umean_int(:) = umean_int(:)/count
    
    Write(*,'(A25,3(2x,g12.5))') 'USF/(1-eps) = ', UMEAN/(one-vol_frac1)
    Write(*,'(A25,3(2x,g12.5))')'UMEAN OF FLUCTUATING FIELD = ', UMEAN_INT
    
    !call tke_bal(ur_tot,umean, mx, my, mz, foffset, nbody, particle_pos,inpos,outpos,energy_bal) 

!!$  tke = 0.d0
!!$  do i=1,mx1
!!$     do k=1,mz
!!$        do j=1,my
!!$           
!!$           tke(1) = tke(1) + (ur_tot(i,j,k,1)-umean_int(1))**2 * particle_pos(i,j,k)
!!$           tke(2) = tke(2) + (ur_tot(i,j,k,2)-umean_int(2))**2 * particle_pos(i,j,k)
!!$           tke(3) = tke(3) + (ur_tot(i,j,k,3)-umean_int(3))**2 * particle_pos(i,j,k)
!!$        end do      !  do j=1,my
!!$     end do         !  do k=1,mz
!!$  end do            !  do i=    
!!$
!!$  tke(:) = tke(:)/dble(count)
  !write(*,'(A25,5D20.12)') 'tke inside fixed bed',tke,sum(tke(:)),sum(tke(:))*0.5
  
  count = 0 
  tke_r = zero 
  do idim=1,ndim
     do k=1,mz
        do j=1,my
           do i=1,mx1
              ur_tot(i,j,k,idim) = (ur_tot(i,j,k,idim)-umean_int(idim)) &
                   & * particle_pos(i,j,k)
              
              arr(i,j,k,idim) = ur_tot(i,j,k,idim)
              tke_r(idim) = tke_r(idim) + arr(i,j,k,idim)**2.d0
              
           end do
        end do
     end do
     call ff3rc(arr(1:mx1,1:my,1:mz,idim), turb_uf(1:my2,1:my,1:mz,idim))
     
  end do

  !turb_uf(:,:,:,1:3) = turb_uf(:,:,:,1:3)*(real(mx1*my*mz,prcn)/real(count_fluid,prcn))
  !count_real = count_fluid!SUM(particle_pos(1:mx1,1:my,1:mz))
  count_real = mx1*my*mz
  write(*,'(A25,5(2x,g12.7))') ' tke real',tke_r/count_real, sum(tke_r(:))/count_real

  tke = 0.d0
  do k=1,mz
     do j=1,my
        do i=1,my/2+1
           if (i==1) then
              tke(1) = tke(1)+dble(turb_uf(i,j,k,1)*conjg(turb_uf(i,j,k,1)))
              tke(2) = tke(2)+dble(turb_uf(i,j,k,2)*conjg(turb_uf(i,j,k,2)))
              tke(3) = tke(3)+dble(turb_uf(i,j,k,3)*conjg(turb_uf(i,j,k,3)))
           else
              tke(1) = tke(1)+dble(turb_uf(i,j,k,1)*conjg(turb_uf(i,j,k,1))) *2.d0
              tke(2) = tke(2)+dble(turb_uf(i,j,k,2)*conjg(turb_uf(i,j,k,2))) *2.d0
              tke(3) = tke(3)+dble(turb_uf(i,j,k,3)*conjg(turb_uf(i,j,k,3))) *2.d0
           end if
        end do
     end do
  end do
  count_comp =1! count_fluid!(mx1*my*mz)**2.d0!1!count_fluid
  write(*,'(A25,5(2x,g12.7))') ' tke complex',tke/dble(count_comp), sum(tke(:))/dble(count_comp)

  !write(*,'(A,4D20.12)') 'tke inside fixed bed(complex)',tke(:),sum(tke(:))
    
  write(*,*) 'difference',sum(tke(:))/(sum(tke_r(:))/(count_real))
!  write(*,*)  size(E_spectra,1)
  E_spectra = 0.d0
  E_uu = 0.d0
  epsf = 0.d0
  do k=1,mz
     do j=1,my
        do i=1,my/2+1
           tmp = -dble(wx(i)*wx(i)+wy(j)*wy(j)+wz(k)*wz(k))
           kmode = sqrt(tmp)
           if (kmode >= kmax) then
              kmode = kmax
           end if
           kindex = int(kmode/dkappa)
           karray(int(kindex)) = kindex*dkappa
           !PRINT*,'kmode = ', kmode, sqrt(tmp)
           tmp = 0.d0
           do idim=1,ndim
              tmp = tmp+dble(turb_uf(i,j,k,idim)*conjg(turb_uf(i,j,k,idim)))
           end do
           tmp1 = dble(turb_uf(i,j,k,1)*conjg(turb_uf(i,j,k,1)))
!           print *,'tmp',tmp
           if (i==1) then
              E_spectra(int(kindex)) = E_spectra(int(kindex)) + tmp
              E_uu(int(kindex)) = E_uu(int(kindex)) + tmp1
              !if (kmode<=22) then
              !   epsf = epsf + tmp1*int(kmode)*int(kmode)
              !end if
           else
              E_spectra(int(kindex)) = E_spectra(int(kindex)) + tmp*2.0
              E_uu(int(kindex)) = E_uu(int(kindex)) + tmp1*2.d0
              !if (kmode<=22) then
              !   epsf = epsf + 2*tmp1*int(kmode)*int(kmode)
              !end if
           end if
        end do
     end do
  end do
!  print *,'after calculate E_spectra'
  print *,'epsf=',epsf
  open(30,file=TRIM(RUN_NAME)//'_E_spectra.dat',form='formatted')
  open(35,file=TRIM(RUN_NAME)//'_E_uu.dat',form='formatted')
  write(30 , *)  'zone t=" total bed  " '
  write(35 , *)  'zone t=" total bed  " '
  ktot = zero 
  do i=0, kindex_max
     !PRINT*,'karrayt = ', karray(i)
     !if(E_spectra(i).gt.zero) then
        write(30,'(I5, 4D20.12)') i, karray(i)*lx/(two*pi), karray(i)*dchar/(two*pi), E_spectra(i)
     
     write(35,'(I5, 4D20.12)') i, dble(i)*eta, dble(i)*dchar, E_uu(i)
  end do
 ! WRITE(*,*)'ktot = ', ktot
  close(30,status="keep")
  deallocate(E_spectra, E_uu,karray)
end SUBROUTINE calc_spec_3d

end MODULE spectra_3d
