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

module init_turb
#include "ibm.h"

	USE precision 
	USE constants
	!use dependent_functions
	use global_data
	use fftw3_interface
	use constants
	use randomno
	use string_funcs
	use nlmainarrays

	type::kappa_type
		real(prcn) :: val
		integer(8) :: num
	end type kappa_type

	!real(prcn), allocatable, dimension(:) :: kx, ky, kz

	type(kappa_type), allocatable :: kappa_tmp(:), kappa(:)
	real(prcn) :: k0, kmin, kmax, dkappa

	real(prcn) :: avr(ndim)
	
	real(prcn) :: eta, l_e, eta_i, l_e_i, eddy_time_i
	real(prcn) :: c_eta, c_l
	character*3  rank_string
	integer :: s_size	

	logical, allocatable :: far_field(:,:,:)

	!real(prcn) :: lx,ly,lz
	real(prcn) :: re_lambda, turb_length, lambda, kmaxeta, u_eta, tau_eta, u_eta_i, tau_eta_i, Re_L, u_prime, epsf_forcing, meaniso_force(ndim), activate_forcing_time, sampling_dissip_time
	integer :: nkappa, mx_iso
	logical :: forced_turbulence

	NAMELIST /turb/ re_lambda, turb_length, nkappa, forced_turbulence, mx_iso, activate_forcing_time, sampling_dissip_time
contains

	subroutine allocation
		implicit none
		real(prcn) :: k2, kmode
		integer :: i, j, k, kindex

		if (I_AM_NODE_ZERO) write (*,"(1a,3d15.7)") "LENGTH OF THE BOX = ", doml

		!	initialize wavenumber vectors
		!allocate(kx(my2), ky(my), kz(my))

		kmin = twopi*1./dble(doml(1))
		k0 = kmin*(1.1) !sqrt(kx(2)**2+ky(2)**2+kz(2)**2)

		if (aliasflag==0) then
			kmax = global_n(1)*kmin/2.
			!kmax = sqrt(2.0)/3.*dble(my)*kmin
		else
			kmax = sqrt(2.0)/3.*dble(my)*kmin
		endif

		dkappa=(kmax-k0)/(nkappa-1)
		allocate(kappa(nkappa), kappa_tmp(nkappa))
		kappa_tmp(:)%num = 0
		kappa(:)%num     = 0



		! Finding number of nodes in each shell
#if !PARALLEL
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
#else
		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
#endif
					if (w2(i,j,k)>small_number) then
						kmode = sqrt(w2(i,j,k))
						if (kmode<k0) then
							kindex=1
						else
							kindex = int((kmode-k0)/dkappa)+2
							if (kindex>nkappa) kindex=nkappa
						endif

#if !PARALLEL
						if (i==1) then
#else
						if (j==1.and.local_o_start(1)==1) then
#endif
							kappa_tmp(kindex)%num = kappa_tmp(kindex)%num + 1
						else
							kappa_tmp(kindex)%num = kappa_tmp(kindex)%num + 2
						endif
					endif
				end do
			end do
		end do

		do i=1, nkappa
			GLOBAL_LONGINT_SUM(kappa_tmp(i)%num,kappa(i)%num,1,comm_cart_2d)
		enddo

		if (I_AM_NODE_ZERO) then
			write (*,*) "SUM OF POINTS IN FOURIER SHELLS = ", sum(kappa(:)%num)
			do i=1,nkappa
				if (kappa(i)%num==0) then
					write (*,"(1a,1i)") "WARNING: THEHE IS A FOURIER SHELL IN WHICH THERE IS NO GRID POINT, IKAPPA = ", i
					write (*,*) "CHANGE THE 'NKAPPA' VALUE"
					STOP
				endif
			enddo
			write (*,*) "ALLOCATION DONE"
		endif

		!INITIALIZING KAPPA%VAL
		do i=1, nkappa
		    kappa(i)%val = kmin+ (i-half)*dkappa
		enddo
	end subroutine allocation

	subroutine calc_initial_turb
		implicit none

		!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
		!IIIII CALCULATING PARAMETERS OF VELOCITY FLUCTUATIONS IIIIIII

		call allocation

		kmin = twopi*1./dble(doml(1))
		k0 = kmin*(1.1) !sqrt(kx(2)**2+ky(2)**2+kz(2)**2)

		if (aliasflag==0) then
			!kmax = global_n(1)*kmin/2.
			kmax = mx*kmin/2.
			!kmax = sqrt(2.0)/3.*dble(my)*kmin
		else
			kmax = sqrt(2.0)/3.*dble(mx)*kmin
		endif


#if 0
		l_e   = turb_length * twopi !doml(1)
		re_l  = re_lambda**2 * 3./20.
		eta   = l_e / ( re_l**0.75 )
		tke   = (re_l * vis / l_e)**2
		epsf  = tke**1.5 / l_e
#endif

		if (.not.(trim(input_type)=='single-phase') .and. mx_iso>0) then
			eta = 2.1 /kmax * dble(mx)/dble(mx_iso)
		!if (trim(input_type)=='single-phase') then
		else
			if (.not.(trim(input_type)=='single-phase')) then
				if (I_AM_NODE_ZERO) write (*,*) "WARNING: MAKE SURE ETA IS SELECTED CORRECTLY"
			endif
			eta = pi/kmax
!eta = dx
		endif
		re_l  = re_lambda**2 * 3./20.
		l_e = eta * ( re_l**0.75 )
		tke   = (re_l * vis / l_e)**2
		epsf  = tke**1.5 / l_e

		!re_lambda = dsqrt(20d0/3d0 * re_l)
		u_prime   = dsqrt(tke * 2d0 / 3d0)
		lambda    = dsqrt(15*vis*u_prime**2/epsf)
		kmaxeta   = kmax*eta
		u_eta	= (epsf*vis)**0.25
		tau_eta	= dsqrt(vis/epsf)

		tke_i  = tke
		epsf_i = epsf
		u_eta_i = u_eta
		tau_eta_i = tau_eta
		eta_i = eta
		l_e_i = l_e
		eddy_time_i = tke_i / epsf_i

		if (I_AM_NODE_ZERO) then
			write (*,*) '-------------------------------'
			write (*,*) 're_lambda = ', re_lambda
			write (*,*) 'Re_L      = ', Re_L
			write (*,*) '-------------------------------'
			write (*,*) 'tke       = ', tke
			write (*,*) 'epsf      = ', epsf
			write (*,*) 'u_prime   = ', u_prime
			write (*,*) 'u_eta     = ', u_eta
			write (*,*) 'tau_eta   = ', tau_eta
			write (*,*) '-------------------------------'
			write (*,*) 'Kmin      = ', kmin
			write (*,*) 'Kmax      = ', kmax
			write (*,*) 'dkappa    = ', dkappa
			write (*,*) 'nkappa    = ', nkappa
			write (*,*) '-------------------------------'
			write (*,*) 'eta       = ', eta
			write (*,*) 'eta/dx    = ', eta/dx
			write (*,*) 'lambda    = ', lambda
			write (*,*) 'L_e       = ', l_e
			write (*,*) 'Kmax*eta  = ', kmaxeta
			write (*,*) '-------------------------------'
		endif
	end subroutine	calc_initial_turb	



	subroutine specinit
!use mypost_process

		use rand_number
		implicit none
		double precision :: kmode
		complex(prcn) :: theta1, theta2
		complex(prcn) :: alpha, beta
		real(prcn)    :: dk

		! local variable
		integer :: i, j, k, m, idim, ierr, unit1, nhbins

		real(prcn) :: k2, umag
		integer :: ikappa,kindex
		real(prcn) :: phi, kxy

		real(prcn) :: tmp
		complex(prcn) :: ctmp
		real(prcn), dimension(:), allocatable :: hist, wt
	
		real(prcn), external :: energy_init, integrate_energy, integrate_epsilon
		real(prcn) :: k_x0, k_mag0
		real(prcn) :: k_x1, k_y1, k_z1
		real(prcn) :: phi_s
		real(prcn) :: u_min,u_max

		character*100 filename

#if PARALLEL
		integer, allocatable :: samplearray(:), which_procs(:)

		complex(prcn), allocatable :: sendarray(:), recvarray(:), plane(:,:,:)

		integer :: mysample, iproc, jproc, ii, jj, jjj, kk, sendid, recvid
#endif

		!--------------------------------------
		if (I_AM_NODE_ZERO) write (*,*) "GENERATING TURBULENT VELOCITY FLUCTUATIONS"

		call calc_initial_turb


		if (I_AM_NODE_ZERO) then
			write (*,*) "FINDING ENERGY SPECTRUM FUNCTION COEFFICIENTS"

			c_eta = 0.4
			c_L   = 6.78

			call xmnewt(c_eta, c_l, tke, epsf_i, l_e, eta, kmax, vis) !for high reynolds numbers

			write (*,"(1a,2d15.7)") "C_eta, C_l = ", c_eta, c_l
			write (*,"(1a,1d17.7)") "ENERGY OF THE SPECTRUM = ", integrate_energy(kmin,kmax,epsf_i, l_e, eta, c_eta, c_l)

			do i=1, nkappa
				if (i==1) then
					k_x0 = 1e-6
					k_x1 = k0
				else
					k_x0 = k0+(i-2) * dkappa
					k_x1 = k0+(i-1) * dkappa
				endif

				! Finding the wave number in each shell (k*) for which the total energy in that shell be euqal to E(k*)*delta_k 
				call xrtnewt(kappa(i)%val,k_x0,k_x1,epsf_i,l_e,eta,c_eta,c_l)
			enddo


			!if (ik_on>0) then
				open(unit=9997,file=trim(run_name)//"_kappa.dat",status="replace",position="append")
				write (9997,*) "zone"
				do i=1,nkappa
					write (9997,"(2D15.7,1I14)") (i-0.5)*dkappa, kappa(i)%val,kappa(i)%num
				enddo
				close(9997)
			!endif
		endif

		BROADCAST_DOUBLE(c_eta,1,node_zero,comm_cart_2d)
		BROADCAST_DOUBLE(c_l,1,node_zero,comm_cart_2d)

		do i=1, nkappa
			BROADCAST_DOUBLE(kappa(i)%val,1,node_zero,comm_cart_2d)
		enddo

		!##################################################
		u(:,:,:,:) = czero

	
		!open(unit=999960,file=trim(run_name)//"_seed.d",status="old",action="read")
		!read (999960,*) seed
		!close(999960)
		!call RLuxGo(389, seed, 0, 0)
		!call RanLux(rand_num1, nsample)
		!call RanLux(rand_num2, nsample)
		!call RanLux(rand_num3, nsample)

		! generating velocity fluctuations
#if !PARALLEL
		do k=1, local_no(3)
			if (I_AM_NODE_ZERO.and.mod(k,10)==0) write (*,*) k,"th PLANE OUT OF ", local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
#else

!		do i=1, local_no(3)
!			do k=1, local_no(2)
!				do j=1, local_no(1)

		do k=1, local_no(2)
			if (I_AM_NODE_ZERO.and.mod(k,10)==0) write (*,*) k,"th PLANE OUT OF ", local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
#endif
					if (w2(i,j,k)<small_number) then
						u(i,j,k,1) = cmplx(0.d0,0.d0)
						u(i,j,k,2) = cmplx(0.d0,0.d0)
						u(i,j,k,3) = cmplx(0.d0,0.d0)
						kmode = 0.d0
						umag = 0.d0
					else
						! The real value of wave number for this node is computed
						k2     = w2(i,j,k)
						kmode  = sqrt(k2)
						k_mag0 = kmode

						! It is determmined with which shell the energy of this node is associated
						if (kmode<k0) then
							kindex=1
						else
							kindex = int((kmode-k0)/dkappa)+2
							if (kindex>nkappa) kindex=nkappa
						endif

						kmode = kappa(kindex)%val

						! The components of the new kappa would be computed to make it consistent with
						! the assumption of the discret energy spectrum on Fourier shells
#if !PARALLEL
						k_x1 = aimag(wx(i))
						k_y1 = aimag(wy(j))
						k_z1 = aimag(wz(k))
#else
						k_x1 = aimag(wy(j))
						k_y1 = aimag(wz(k))
						k_z1 = aimag(wx(i))
#endif

						phi_s = kmode/k_mag0
						k_x1  = k_x1*phi_s
						k_y1  = k_y1*phi_s
						k_z1  = k_z1*phi_s

						kxy   = dsqrt(k_x1**2+k_y1**2)
						!------------------------

						call gen_rand

						theta1 = cmplx(0,twopi*rand1)
						theta2 = cmplx(0,twopi*rand2)
						phi    = twopi*rand3

						umag = energy_init(kmode, epsf_i, l_e_i, eta_i, c_eta, c_l)

						! alpha and beta are computed according to the method presented
						! E(k)*dkappa produces the energy of each shell. Then it is distributed
						! uniformly among nodes related to the shell. Acoordingly, the manitude of
						! each velocity fluctuation is computed according to the amount of energy
						! it has.

						dk=dkappa
						if (kindex==1) dk=k0

						alpha = sqrt(2 * umag * dk / kappa(kindex)%num) * exp(theta1) * cos(phi)
						beta  = sqrt(2 * umag * dk / kappa(kindex)%num) * exp(theta2) * sin(phi)


						if (abs(kxy)<=small_number) then
							u(i,j,k,1) = alpha
							u(i,j,k,2) = beta
							u(i,j,k,3) = czero
						else
							u(i,j,k,1) = alpha*kmode*k_y1+beta*k_x1*k_z1
							u(i,j,k,1) = u(i,j,k,1)/(kmode*kxy)

							u(i,j,k,2) = beta*k_y1*k_z1-alpha*kmode*k_x1
							u(i,j,k,2) = u(i,j,k,2)/(kmode*kxy)

							u(i,j,k,3) = -1d0*beta*kxy/kmode
						endif
					endif
				enddo
			enddo
		enddo

#if !PARALLEL
		!for i=1,j=1,k<>1 line
		i=1
		j=1
		do k=mz/2+1,mz
			u(i,j,k,1) = conjg(u(i,j,mz+2-k,1))
			u(i,j,k,2) = conjg(u(i,j,mz+2-k,2))
			u(i,j,k,3) = conjg(u(i,j,mz+2-k,3))
		enddo

		!for i=1, j<>1, k=1 line
		i=1
		k=1
		do j=2,my/2
			u(i,my+2-j,k,1) = dcmplx(dble(u(i,j,k,1)),(-1.0)*aimag(u(i,j,k,1)))
			u(i,my+2-j,k,2) = dcmplx(dble(u(i,j,k,2)),(-1.0)*aimag(u(i,j,k,2)))
			u(i,my+2-j,k,3) = dcmplx(dble(u(i,j,k,3)),(-1.0)*aimag(u(i,j,k,3)))
		enddo

		!for i=1, J<>1, k<>1 
		 i=1
		 do k=2, my/2
		    do j=2, my/2
		       u(i,my+2-j,mz+2-k,:) = conjg(u(i,j,k,:))
		       u(i,my+2-j,k,:)      = conjg(u(i,j,mz+2-k,:))
		    end do
		 enddo
#else



		allocate(samplearray(0:nproc-1))

		mysample = local_no(3)*local_no(2)*ndim
		call mpi_allgather(mysample, 1, mpi_int, samplearray, 1, mpi_int, comm_cart_2d, ierr)

		allocate(which_procs(nprocz))
		which_procs = -1

		jproc=0
		do iproc=0, nproc-1
			if (xstarts(iproc)==1 .and. xlocals(iproc)>0  ) then
				jproc = jproc+1

				if (jproc>nprocz) then
					if (I_AM_NODE_ZERO) write (*,*) "ERROR IN FORMING WHICH_PROCS"
					PARALLEL_FINISH()
					stop
				endif
				which_procs(jproc) = iproc
			endif
		enddo

		if (local_o_start(1)==1 .and. local_no(1)>0) then
			allocate(plane(global_n(3),global_n(2),ndim))
			plane = czero

			allocate(sendarray(mysample))
			sendarray = czero

			m = 0
			j  = 1
			do idim=1, ndim
				do k=1, local_no(2)
					do i=1, local_no(3)
						ii = local_o_start(3)+i-1
						kk = local_o_start(2)+k-1

						plane(ii,kk,idim) = u(i,j,k,idim)

						m = m+1
						sendarray(m) = u(i,j,k,idim)
					enddo
				enddo
			enddo

			if (nprocz>1) then
				do iproc=1, nprocz
					sendid = which_procs(iproc)

					if (myid==sendid) then
						do jproc=1, nprocz
							recvid = which_procs(jproc)
							if (sendid/=recvid) then

								call mpi_send(sendarray, mysample, mpi_double_complex, recvid, sendid, comm_cart_2d, ierr)
							endif
						enddo

					else
						allocate(recvarray(samplearray(sendid)))
						recvarray = czero
					
						recvid = myid
						if (sendid/=recvid) then

							call mpi_recv(recvarray, samplearray(sendid), mpi_double_complex, sendid, sendid, comm_cart_2d, status, ierr)
						endif

						m = 0
						do idim=1, ndim
							do k=ystarts(sendid), yends(sendid)
								do i=zstarts(sendid), zends(sendid)
									m = m+1

									plane(i,k,idim) = recvarray(m)
								enddo
							enddo
						enddo

						deallocate(recvarray)
					endif
				enddo
			endif

			deallocate(sendarray)

			!for i=1,j=1,k<>1 line
			if (local_o_start(1)==1 .and. local_o_start(2)==1) then
				i=1
				j=1
				do k=mz/2+1,mz
					u(k,i,j,1) = conjg(plane(mz+2-k,j,1))
					u(k,i,j,2) = conjg(plane(mz+2-k,j,2))
					u(k,i,j,3) = conjg(plane(mz+2-k,j,3))
				enddo
			endif


			!for i=1, j<>1, k=1 line
			if (local_o_start(1)==1 .and. local_o_start(3)==1) then
				i=1
				k=1
				do j=2,my/2
					jj = my+2-j

					if (local_o_start(2)<=jj .and. jj<=local_o_end(2)) then
						jjj = jj - local_o_start(2) + 1
				
						u(k,i,jjj,1) = dcmplx(dble(plane(k,j,1)),(-1.0)*aimag(plane(k,j,1)))
						u(k,i,jjj,2) = dcmplx(dble(plane(k,j,2)),(-1.0)*aimag(plane(k,j,2)))
						u(k,i,jjj,3) = dcmplx(dble(plane(k,j,3)),(-1.0)*aimag(plane(k,j,3)))
					endif
				enddo
			endif


			!for i=1, J<>1, k<>1 
			if (local_o_start(1)==1) then
				i=1
				do k=2, my/2
					do j=2, my/2
						jj = my+2-j

						if (local_o_start(2)<=jj .and. jj<=local_o_end(2)) then
							jjj = jj - local_o_start(2)+1 

							u(mz+2-k,i,jjj,:) = conjg(plane(k,j,:))
							u(k,i,jjj,:)      = conjg(plane(mz+2-k,j,:))
						endif
					 enddo
				 enddo
			endif
		endif
#endif

!^^^^ existance of this part does not have any effect on the answers
!		u(1,1,my/2+1,:) = czero
!		u(1,ny/2+1,1,:) = czero
!		u(1,ny/2+1,:,1) = czero
!		u(1,ny/2+1,:,2) = czero
!		u(1,ny/2+1,:,3) = czero
!
!		u(1,:,nz/2+1,1) = czero
!		u(1,:,nz/2+1,2) = czero
!		u(1,:,nz/2+1,3) = czero
!
!		u(nx/2+1,:,:,:) = czero     ! i=nxturbby2p1 plane is set to zero
!		if they are added in then, some of the energy would be lost
!-------------------------------------------------

		call calc_velreal(u, umean, ubcp)

		call compute_average(0)
		call statistics(0)

		if (I_AM_NODE_ZERO) then
			write (*,*) "SPECINIT.F90 IS DONE..."
			call screen_separator(80,'I')
		endif

#if 0
open (1,file="u.bin", status="replace", action="write", form="formatted")
open (2,file="v.bin", status="replace", action="write", form="formatted")
open (3,file="w.bin", status="replace", action="write ", form="formatted")

do k=1, global_n(3)
	do j=1, global_n(2)
		do i=1, global_n(1)
			write (1,*) ubcp(i,j,k,1)
			write (2,*) ubcp(i,j,k,2)
			write (3,*) ubcp(i,j,k,3)
		enddo
	enddo
enddo
close(1)
close(2)
close(3)


open (1,file="u.bin", status="old", action="read", form="formatted")
open (2,file="v.bin", status="old", action="read", form="formatted")
open (3,file="w.bin", status="old", action="read", form="formatted")

do k=1, global_n(3)
	do j=1, global_n(2)
		do i=1, global_n(1)
			read (1,*) ubcp(i,j,k,1)
			read (2,*) ubcp(i,j,k,2)
			read (3,*) ubcp(i,j,k,3)
		enddo
	enddo
enddo
close(1)
close(2)
close(3)

do idim=1, ndim
	call fftwr2c(ubcp(1:local_ni(1),1:local_ni(2),1:local_ni(1),idim), u(:,:,:,idim))
enddo
write (*,*) "READING VELOCITY FIELD FROM SAVE FILES"
#endif

	end subroutine specinit

	subroutine compute_average(call_flag)
		implicit none

		integer :: call_flag
		!real(prcn), dimension(:,:,:,:), pointer :: u_r

		real(prcn) :: tmp(3)
		integer :: i, j, k, idim
		!integer :: imax

		!u_r  => ubcp
		!imax =  nx

		avr = zero
		tmp = zero
		do idim=1, ndim
			do k=1, local_ni(3)
				do j=1, local_ni(2)
					do i=1, local_ni(1)
						tmp(idim) = tmp(idim) + ubcp(i,j,k,idim)
					enddo
				enddo
			enddo
		enddo

		GLOBAL_DOUBLE_SUM(tmp,avr,3,comm_cart_2d)

		avr(:) = avr(:) / dble(global_n(1))/dble(global_n(2))/dble(global_n(3)) 
		if (I_AM_NODE_ZERO) write (*,"(1a,3d15.7)") "UAVERAGE = ", avr(:)
	end subroutine compute_average

	subroutine statistics(call_flag)
		use general_funcs, only : instant_file_opener
		implicit none
		integer :: call_flag

		real(prcn) :: kmode, k2, u_prime, k_x0, k_x1, dk
		real(prcn), external :: energy_init, integrate_energy, integrate_epsilon

		real(prcn), dimension(:), allocatable :: engmag, engmag_tmp, epsfmag, epsfmag_tmp
		real(prcn), dimension(:), allocatable :: hist, wt

		!real(prcn),    dimension(:,:,:,:), pointer :: u_r
		!complex(prcn), dimension(:,:,:,:), pointer :: u_f

		real(prcn) :: krange, energyf_pre(ndim)
		real(prcn) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8

		integer :: i, j, k, l, m, idim, ierr, unit1, nhbins, nsample
		integer :: imax, jmax, kkmax, ikrange
		integer	:: ikappa,kindex

		real(prcn) :: l_e, Re_L
		real(prcn) :: umag, phi

		real(prcn) :: var(ndim), skew(ndim), kurt(ndim), turb(ndim), bij(ndim), vectmp(ndim)
		real(prcn) :: u_min, u_max, tmp, tmp1, tmp2, tmp3, mean_energy
		complex(prcn) :: ctmp
		integer :: unitno
		logical :: filexist
		character*100 filename1
		character*7  step_string
		integer :: filenum
		!-------------------------------------------------------------

		!u_f => u
		!u_r => ubcp

		!call check_divergence

		if (.not.allocated(kappa)) call allocation

		tke = 0d0
		epsf = 0d0
		energyf_pre = 0.d0

		allocate(engmag(nkappa),engmag_tmp(nkappa))
		allocate(epsfmag(nkappa),epsfmag_tmp(nkappa))
		engmag = zero
		engmag_tmp = zero

		epsfmag = zero
		epsfmag_tmp = zero

#if !PARALLEL
		do k=1, local_no(3)
			if (mod(k,10)==0) write (*,*) k,"th PLANE OUT OF ", my
			do j=1, local_no(2)
				do i=1, local_no(1)
#else
		do k=1, local_no(2)
			if (I_AM_NODE_ZERO .and. mod(k,10)==0) write (*,*) k,"th PLANE OUT OF ", local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
#endif
					k2    = w2(i,j,k)
					kmode = dsqrt(k2)
					if (kmode<k0) then
						kindex=1
					else
						kindex = int((kmode-k0)/dkappa)+2
						if (kindex>nkappa) kindex=nkappa
					endif

					kmode=kappa(kindex)%val

					tmp1 = dble(u(i,j,k,1) * conjg(u(i,j,k,1)))
					tmp2 = dble(u(i,j,k,2) * conjg(u(i,j,k,2)))
					tmp3 = dble(u(i,j,k,3) * conjg(u(i,j,k,3)))

#if !PARALLEL
					if ( i==1 ) then
#else
					if ( j==1 .and. local_o_start(1)==1 ) then
#endif
						tmp1 = tmp1/2
						tmp2 = tmp2/2
						tmp3 = tmp3/2
					endif
					energyf_pre(1) = energyf_pre(1)+tmp1
					energyf_pre(2) = energyf_pre(2)+tmp2
					energyf_pre(3) = energyf_pre(3)+tmp3

					tmp = tmp1+tmp2+tmp3

					tke = tke + tmp
					epsf = epsf + 2 * vis * k2 * tmp

					engmag_tmp(kindex) = engmag_tmp(kindex) + tmp
					epsfmag_tmp(kindex)= epsfmag_tmp(kindex)+ 2*vis*k2*tmp
				enddo
			enddo
		enddo

		if (call_flag==0) then
			tmp1 = tke
			GLOBAL_DOUBLE_SUM(tmp1,tke,1,comm_cart_2d)

			tmp1 = epsf
			GLOBAL_DOUBLE_SUM(tmp1,epsf,1,comm_cart_2d)
		endif

		vectmp(:) = energyf_pre(:)
		GLOBAL_DOUBLE_SUM(vectmp,energyf_pre,3,comm_cart_2d)

		GLOBAL_DOUBLE_SUM(engmag_tmp,engmag,nkappa,comm_cart_2d)
		GLOBAL_DOUBLE_SUM(epsfmag_tmp,epsfmag,nkappa,comm_cart_2d)

10		if (I_AM_NODE_ZERO) then
			write (*,"(1a,2d15.7)") "ENERGY SUM = ", sum(engmag(:)), tke
			write (*,"(1a,2d15.7)") "DISSIP SUM = ", sum(epsfmag(:)), epsf

			eta  = (vis**3 / epsf)**(0.25d0)
			l_e  = tke**1.5 / epsf
			re_l = l_e * dsqrt(tke) / vis
			re_lambda = dsqrt(20d0/3d0 * re_l)
			u_prime   = dsqrt(tke * 2d0 / 3d0)
			lambda    = dsqrt(15*vis*u_prime**2/epsf)
			kmaxeta   = kmax*eta
			u_eta	= (epsf*vis)**0.25
			tau_eta	= dsqrt(vis/epsf)
		
20			write (*,*) "###### OUTPUT ACCORDING TO SINGLE PHASE FLOW #####"
			write (*,*) '-------------------------------'
			write (*,*) 'C_eta     = ', c_eta
			write (*,*) 'C_L       = ', c_l
			write (*,*) '-------------------------------'
			!write (*,*) 'Re_p      = ', re
			write (*,*) 're_lambda = ', re_lambda
			write (*,*) 'Re_L      = ', Re_L
			write (*,*) '-------------------------------'
			write (*,*) 'tke       = ', tke
			write (*,*) 'epsf      = ', epsf
			write (*,*) 'vis       = ', vis
			write (*,*) 'u_prime   = ', u_prime
			write (*,*) 'u_eta     = ', u_eta
			write (*,*) 'tau_eta   = ', tau_eta
			write (*,*) '-------------------------------'
			write (*,*) 'Kmin      = ', kmin
			write (*,*) 'Kmax      = ', kmax
			write (*,*) 'dkappa    = ', dkappa
			write (*,*) 'nkappa    = ', nkappa
			write (*,*) '-------------------------------'
			write (*,*) 'eta       = ', eta
			write (*,*) 'dx/eta    = ', dx / eta
			write (*,*) 'eta/L_e   = ', eta / l_e
			write (*,*) 'eta/L     = ', eta / doml(2)
			write (*,*) 'lambda    = ', lambda
			write (*,*) 'lambda/dx = ', lambda / dx
			write (*,*) 'lambda/L_e= ', lambda / l_e
			write (*,*) 'lambda/L  = ', lambda / doml(2)
			write (*,*) 'L_e       = ', l_e
			write (*,*) 'L_e/dx    = ', l_e/dx
			write (*,*) 'L_e/eta   = ', l_e/eta
			write (*,*) "Box Length= ", doml(2)
			write (*,*) 'L_e/Ly    = ', l_e / doml(2)
			write (*,*) 'Kmax*eta  = ', kmaxeta
			write (*,*) '-------------------------------'
			write (*,*) 'eta/eta_i         = ', eta/eta_i
			write (*,*) 'l_e/l_e_i         = ', l_e/l_e_i
			write (*,*) 'tke/tke_i         = ', tke/tke_i
			write (*,*) 'eps/eps_i         = ', epsf/epsf_i
			write (*,*) 'u_eta/u_eta_i     = ', u_eta/u_eta_i
			write (*,*) 'tau_eta/tau_eta_i = ', tau_eta/tau_eta_i

			write (*,*) '-------------------------------'
			if (kmaxeta<1.5) write (*,*) "WARNING: KAMX*ETA<1.5; TE CRITERION IS NOT SATISFIED"

			filename1 = trim(run_name)//"_statistics"
			filenum = 1
			call instant_file_opener(filename1,filenum,.true.)

			write (filenum,"(21D15.7)") t/t_conv, re_lambda, re_l, tke, u_prime, epsf, eta, eta/l_e, eta/doml(2), &
						&  lambda, lambda/l_e, lambda/doml(2), l_e, l_e/doml(2), u_eta, &
						&  tau_eta, kmaxeta, c_eta, c_l
			close (filenum)
			call screen_separator(80,'S')


			filename1 = trim(run_name)//"_spectra"
			filenum = 1
			call instant_file_opener(filename1,filenum,.true.)

			if (call_flag==0) then
				write (filenum,*) "zone"
				!integration of energy spectrum function: k0-kmax
				krange = kappa(1)%val
				dk = 0.01
				ikrange = int((kmax-krange)/dk)

				do i=1,ikrange
					sum3=integrate_energy(krange+(i-1)*dk,krange+i*dk,epsf,l_e,eta,c_eta,c_l)

					kmode = krange +i*dk
					umag  = energy_init(kmode,epsf,l_e,eta,c_eta,c_l)
					write(filenum,'(1I,9D15.7)') i,kmode, kmode*eta, kmode*l_e,	&
							& umag, umag/(eta*u_eta*u_eta),umag/(tke*l_e),		&
							& 2.0*umag*vis*kmode**2, 2.0*umag*vis*kmode**2/epsf,sum3/tke
				enddo
			endif

			!integration of discrete energy spectrum function: k0-kmax
			write (filenum,*) "zone"
			sum1=0d0
			do i=1,nkappa
				kmode = kappa(i)%val
				tmp	 = engmag(i)
				phi  = epsfmag(i)
		
				sum1 = sum1 + tmp
		
				if (i==1) then
					tmp = tmp / k0
					phi = phi / k0
				else
					tmp = tmp / dkappa
					phi = phi / dkappa
				endif

				write(filenum,'(1I,9D15.7)') i, kmode, kmode*eta, kmode*l_e,     &
						& tmp, tmp/(eta*u_eta*u_eta), tmp/(tke*l_e),			&
						& phi, phi/epsf, engmag(i)/tke
			enddo
			close(filenum)

			if (call_flag==0) then
				!integration of energy spectrum function: 0-k0
				sum2 = integrate_energy(1d-6,k0,epsf,l_e,eta,c_eta,c_l)

				!integration of energy spectrum function: k0-kmax
				sum3 = integrate_energy(k0,kmax,epsf,l_e,eta,c_eta,c_l)

				!integration of energy spectrum function: kmax-infinity
				sum4 = integrate_energy(kmax,3*kmax,epsf,l_e,eta,c_eta,c_l)

				!integration for epsilon: 0-k0
				sum5 = integrate_epsilon(1d-6,k0,vis,epsf,l_e,eta,c_eta,c_l)

				!integration for epsilon: k0-kmax
				sum6 = integrate_epsilon(k0,kmax,vis,epsf,l_e,eta,c_eta,c_l)

				!integration for epsilon: kmax-infinity
				sum7 = integrate_epsilon(kmax,3*kmax,vis,epsf,l_e,eta,c_eta,c_l)
	
				write (*,"(A,D)") 'TKE from Fourier space     (1) = ', tke
				write (*,"(A,D)") 'TKE from E.S.F,    0->k0   (2) = ', sum2
				write (*,"(A,D)") 'TKE from E.S.F,   k0->kmax (3) = ', sum3
				write (*,"(A,D)") 'TKE from E.S.F, kmax->inf. (4) = ', sum4
				write (*,"(A,D)") 'TKE from E.S.F,    0->kmax (5) = ', sum2 + sum3
				write (*,"(A,D)") 'TKE INITIAL                    = ', tke_i
				call screen_separator(40,'.')
	
				write (*,"(A,D)") '(1)/(5)                        = ', tke / (sum2 + sum3)
				write (*,"(A,D)") '(1)/TKE_I                      = ', tke / tke_i
				write (*,"(A,D)") '(5)/TKE_I                      = ', (sum2 + sum3) / tke_i

				call screen_separator(50,'*')

				write (*,"(A,D)") 'EPS from Fourier space     (6) = ', epsf
				write (*,"(A,D)") 'EPS from E.S.F,    0->k0   (7) = ', sum5
				write (*,"(A,D)") 'EPS from E.S.F,   k0->kmax (8) = ', sum6
				write (*,"(A,D)") 'EPS from E.S.F, kmax->inf. (9) = ', sum7
				write (*,"(A,D)") 'EPS from E.S.F,    0->kmax (10)= ', sum5 + sum6
				write (*,"(A,D)") 'EPS INITIAL                    = ', epsf_i
				call screen_separator(40,'.')

				write (*,"(A,D)") '(6) /(10)                      = ', epsf / (sum5 + sum6)
				write (*,"(A,D)") '(6) /EPS_I                     = ', epsf / epsf_i
				write (*,"(A,D)") '(10)/EPS_I                     = ', (sum5+sum6) / epsf_i
				call screen_separator(80,'I')
				!-------------------------------------------------
				! checking phase space velocity in fluid space
				!-------------------------------------------------
			endif
		endif
	
#if 0
		if (pdf_on) then
			unit1 = 9998
			open (unitno,file=trim(run_name)//"_pdf.dat", status="replace")
			write (unit1,*) "variables=V,F(V)"
			write (unit1,*) 'Zone t="',0,'"'

			j=50
			do i=-j,j
				phi = 5d0/j * i
				tmp = 1d0/dsqrt(2*pi)*dexp(-phi**2/2)
				write (unit1,"(2D15.7)") phi, tmp
			enddo
		endif
#endif


		if (call_flag==0) then
			call gauss3D(avr(:), var(:), skew(:), kurt(:))	
			do idim = 1, ndim
				u_max = maxval(ubcp(:,:,:,idim))
				u_min = minval(ubcp(:,:,:,idim))

				tmp1 = u_max
				GLOBAL_DOUBLE_MAX(tmp1,u_max,1,comm_cart_2d)

				tmp1 = u_min
				GLOBAL_DOUBLE_MIN(tmp1,u_min,1,comm_cart_2d)

				turb(idim) = var(idim) / 2
				if (I_AM_NODE_ZERO) then
					write (*,*) "For direction # ", idim
					write (*,*) "u_min/r.m.s. = " , u_min /dsqrt(var(idim))
					write (*,*) "u_max/r.m.s. = " , u_max /dsqrt(var(idim))
					write (*,*) "Average     = " , avr(idim)
					write (*,*) "Variance    = " , var(idim)
					!write (*,*) "Skewness    = " , skew
					!write (*,*) "Flatness    = " , kurt
					write (*,'(A,2D15.5)') "0.5<u*u> (physical, Fourier space)= " , turb(idim), energyf_pre(idim)!,turb()-energyf_pre(1) 
					write (*,*) "..............................."
					write (*,*) "-------------------------------"
				endif

#if 0
				if (pdf_on) then
					if (.not.allocated(wt)) then
						nsample = my * my * my
						allocate(wt(nsample))
						do i=1, nsample
							wt(i)=1d0/nsample
						enddo
						nhbins=50
						allocate(hist(nhbins))
					endif
			
					CALL histogram1(ubcp(:,:,:,idim),wt,nsample,nhbins,u_min,u_max,hist(1:nhbins))
					CALL plothist(hist(1:nhbins),u_min,u_max,nhbins,unit1,real(idim,prcn),1.d0)
					if (idim==3) close(unit1)
				endif
#endif
			enddo
			if (I_AM_NODE_ZERO) write (*,"(A,1D)") 'TKE from inverse of the Fourier space = ', sum(turb(1:3))
		endif	

		if (allocated(wt)) then
			deallocate(wt)
			deallocate(hist)
		endif
		deallocate(engmag,epsfmag)
		deallocate(engmag_tmp,epsfmag_tmp)
		!-------------------------------------------------
	 
#if 0
		write (*,*) "WRITING THE VELOCITY FLUCTUATIONS INTO OUTPUT FILE"
		open(unit=9994,file=trim(run_name)//"_velocity.dat",status="replace",form="unformatted")
		write (9994) re_lambda, Re_L, tke, epsf, vis, u_eta, tau_eta, kmin, kmax, eta, lambda, l_e, kmaxeta

		do k=1, mz
			do j=1, my
				do i=1,nx
					write (9994) i, j, k, u_r(i,j,k,:)
				enddo
			enddo
		enddo
		close (9994)
#endif
	end subroutine statistics

	subroutine gauss3D(ave, var, skew, kurt)
		implicit none
		double precision, intent(inout) :: ave(:), var(:), skew(:), kurt(:)

		!------------------------
		! local variables

		integer :: i,j,k,ierr
		double precision :: tmp(ndim), tkein(ndim)

		tmp(:) = zero
		avr(:) = zero
		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					tmp(:) = tmp(:) + ubcp(i,j,k,:)
				enddo
			enddo
		enddo

		GLOBAL_DOUBLE_SUM(tmp,ave,3,comm_cart_2d)
		ave(:) = ave(:)/dble(global_n(1))/dble(global_n(2))/dble(global_n(3))
	  
		var(:) = 0.d0
		tkein(:) = zero
		skew(:) = 0.d0
		kurt(:) = 0.d0
		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					var(:) = var(:) + (ubcp(i,j,k,:)-ave(:))*(ubcp(i,j,k,:)-ave(:))  ! var
					!tkein(:) = tkein(:) + ubcp(i,j,k,:)*ubcp(i,j,k,:)              
					skew(:) = skew(:) + (ubcp(i,j,k,:)-ave(:))**3                ! skew
					kurt(:) = kurt(:) + (ubcp(i,j,k,:)-ave(:))**4                ! kurt
				enddo
			enddo
		enddo

		tmp(:) = var(:)
		var(:) = zero
		GLOBAL_DOUBLE_SUM(tmp,var,3,comm_cart_2d)

		tmp(:) = skew(:)
		skew(:) = zero
		GLOBAL_DOUBLE_SUM(tmp,skew,3,comm_cart_2d)
		
		tmp(:) = kurt(:)
		kurt(:) = zero
		GLOBAL_DOUBLE_SUM(tmp,kurt,3,comm_cart_2d)

		var  = var / dble(global_n(1))/dble(global_n(2))/dble(global_n(3))
		!tkein = tmp_r(2)/dble(nx*ny*nz_tot)*0.5
		skew = skew / dble(global_n(1))/dble(global_n(2))/dble(global_n(3)) / (var**1.5) !/ var**1.5
		kurt = kurt / dble(global_n(1))/dble(global_n(2))/dble(global_n(3)) / (var**2) !/ var**2
	end subroutine gauss3D


	subroutine linear_forcing
		use general_funcs, only : instant_file_opener
		use bcsetarrays
		implicit none

		integer :: i, j, k, idim
		real(prcn) :: force_tmp(ndim)
		character*100 filename
		integer :: unitnum
		if (trim(input_type).eq."single-phase") fr(:,:,:,:) = zero

		if (I_AM_NODE_ZERO) then
			call screen_separator(80,'F')
			write (*,*) "FORCING TERM = ", epsf_forcing/(2*tke)
		endif

		do idim=1, ndim
			do k=1, local_ni(3)
				do j=1, local_ni(2)
					do i=1, local_ni(1)
						if (fluid_atijk(i,j,k)) then
							fr(i,j,k,idim) = (ubcp(i,j,k,idim)-ufmean(idim)) * epsf_forcing/(2*tke)
						endif
					enddo
				enddo
			enddo
		enddo


		meaniso_force(:) = zero
		force_tmp(:) = zero
		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					if (fluid_atijk(i,j,k)) then
						force_tmp(:) = force_tmp(:) + fr(i,j,k,:)
					endif
				enddo
			enddo
		enddo
		GLOBAL_DOUBLE_SUM(force_tmp,meaniso_force,3,comm_cart_2d)

		meaniso_force(:) = meaniso_force(:) / count_fluid

		if (I_AM_NODE_ZERO) then
			write (*,"(1a,3d15.7)") "MEAN ISOTROPIC FORCING = ", meaniso_force(:)

			filename=trim(run_name)//"_forcing_term"
			unitnum=1
			call instant_file_opener(filename,unitnum,.true.)

			write (unitnum,"(2d15.7)") t/t_conv, epsf_forcing/(2*tke)
			close(unitnum)
			call screen_separator(80,'F')
		endif

	end subroutine linear_forcing


	subroutine calc_velreal(velc, velmean, velr)
		use parallel
		implicit none 

		complex(prcn), intent(inout) :: velc(:,:,:,:)
		real(prcn), intent(in) :: velmean(1:ndim)

		real(prcn), intent(inout) :: velr(-1:local_ni(1)+2,-1:local_ni(2)+2,-1:local_ni(3)+2,1:ndim)
		integer :: idim

		do idim=1, ndim
			call fftwc2r(velc(:,:,:,idim), &
						&	 velr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim))

			velr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim) = &
					&	velr(1:local_ni(1),1:local_ni(2),1:local_ni(3),idim) + velmean(idim)
		enddo
		call communicate_velocity(velr(-1:local_ni(1)+2,-1:local_ni(2)+2,-1:local_ni(3)+2,1:ndim))
	end subroutine calc_velreal

	subroutine check_divergence
		IMPLICIT NONE
		INTEGER :: i, j, k, idim
		COMPLEX(prcn) :: tmpc, wtmp
		REAL(prcn) :: tmp, divmax
    
		divmax = zero
#if !PARALLEL
		do k = 1, local_no(3)
			do j = 1, local_no(2)
				do i = 1, local_no(1)
					tmpc = czero
					do idim=1, ndim

						if (idim == 1) then
							wtmp = wx(i)
						elseif (idim == 2) then
							wtmp = wy(j)
						elseif (idim == 3) then
							wtmp = wz(k)
						endif
#else
		do k = 1, local_no(2)
			do j = 1, local_no(1)
				do i = 1, local_no(3)
					tmpc = czero
					do idim=1, ndim

						if (idim == 1) then
							wtmp = wy(j)
						elseif (idim == 2) then
							wtmp = wz(k)
						elseif (idim == 3) then
							wtmp = wx(i)
						endif
#endif
						tmpc = tmpc + u(i,j,k,idim)*wtmp
					enddo
					tmp = DSQRT(dble(tmpc*conjg(tmpc)))

					if (tmp.gt.divmax) divmax = tmp
				enddo 
			enddo
		enddo
		GLOBAL_DOUBLE_MAX(divmax, tmp, 1, comm_cart_2d)

		if (I_AM_NODE_ZERO) write (*,'(A25,g12.5,3i)') 'MAX DIVERGENCE = ', tmp
	end subroutine check_divergence


end module init_turb

SUBROUTINE histogram1(u,wt,n,nddu,bldd,brdd,hist)
  !nddu: number of bins for pdf formation
  !bldd: lefh hand side limit .... output
  !brdd: right side limie.... otuput
  USE precision
  USE constants
  USE randomno
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n , nddu
  REAL(prcn), INTENT(in), DIMENSION(n) :: u, wt
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: u_t
  
  REAL(prcn),  INTENT(out) :: bldd, brdd 

  
  REAL(prcn),  INTENT(inout), DIMENSION(nddu) ::  hist
    REAL(prcn) ::  xdiff
    
    REAL(prcn) :: vlmt, ave, adev, sdev, var, skew, curt
    
    INTEGER :: i, ibin
    bldd= 1.e25
    brdd=-1.e25

    ALLOCATE(u_t(n))
    CALL moment1(4, n, n, wt, u, ave,adev,sdev,var,skew&
         &,curt)
!    PRINT*,'ave,var,sdev, skew, curt=', ave, var,sdev, skew,curt
!    WRITE(*,*)'number of bins in hist..',nddu

  DO i=1,n
     u_t(i)=(u(i)-ave)/sdev
  ENDDO

  CALL moment1(4, n, n, wt, u_t, ave,adev,sdev,var,skew&
       &,curt)
!!$  PRINT*,'ave,var, sdev,skew, curt=', ave, var,sdev, skew,curt

    DO i=1,nddu
       hist(i)=0.0  
    ENDDO
    bldd = MIN(MINVAL(u_t(:)), bldd)
    brdd = MAX(MAXVAL(u_t(:)), brdd)
    xdiff = (brdd - bldd)/float(nddu-1)

    DO i=1,n
       ibin = (u_t(i) - bldd)/xdiff + 1
       hist(ibin)=hist(ibin) + wt(i)/xdiff
    ENDDO

    DEALLOCATE(u_t)
END SUBROUTINE histogram1


SUBROUTINE plothist(hist,lb,ub,nbins,iunit,t,tref)
  USE precision
  USE constants
  IMPLICIT NONE
  REAL(prcn), INTENT(in), DIMENSION(nbins) ::hist
  
  INTEGER, INTENT(in) ::iunit,nbins
  REAL(prcn) ::  sum_x,lb,ub,dx, t, tref, tmp, tmp2
  INTEGER(prcn) :: i
  
  sum_x = lb
  dx= (ub-lb)/(float(nbins)-1)
  tmp = one/sqrt(twopi)
  WRITE(iunit,*)'Zone t="',int(t/tref),'"'
  DO i=1,nbins
     tmp2  = sum_x+dx/2
     WRITE(iunit,*)tmp2,hist(i)!, tmp*exp(-(tmp2*tmp2)/two)
     sum_x = sum_x + dx
  ENDDO
  
  RETURN
END SUBROUTINE plothist


real(prcn) function energy_init(kmode,epsf,l_e,eta,c_eta,c_l)
	use precision

	real(prcn), intent(in) :: kmode, epsf, l_e, eta, c_eta, c_l
	real(prcn), parameter  :: C = 1.5, beta = 5.2, p0 = 2.0

!	c_eta = 0.4,c_L = 6.78,

	energy_init = C * epsf**(2./3.) * kmode**(-5./3.)
	energy_init = energy_init*(kmode * l_e/sqrt((kmode*l_e)**2 +&
								& c_L))**(5./3. + p0)
	energy_init = energy_init*exp(-1.d0*beta*(((kmode*eta)**4 + c_eta**4)&
								&**0.25 - c_eta))
	return   
end function energy_init

real(prcn) function integrate_energy(kmin,kmax,eps,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) :: kmin,kmax,eps,l_e,eta,c_eta,c_l
	
	real(prcn), external :: energy_init
	real(prcn) :: dk, umag0, umag, kmode
	integer :: i,ikrange

	dk = 0.001
	ikrange = int((kmax-kmin)/dk)

	integrate_energy  = 0d0
	umag0 = energy_init(kmin,eps,l_e,eta,c_eta,c_l)
	do i=1,ikrange
		kmode            = kmin +i*dk
		umag             = energy_init(kmode,eps,l_e,eta,c_eta,c_l)
		integrate_energy = integrate_energy + (umag+umag0)*dk/2
		umag0            = umag
    end do
end function integrate_energy

real(prcn) function integrate_epsilon(kmin,kmax,vis,epsf,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) :: kmin,kmax,vis,epsf,l_e,eta,c_eta,c_l
	
	real(prcn), external :: energy_init
	real(prcn)	:: dk, umag0, umag, kmode
	integer :: i,ikrange

	dk = 0.001
	ikrange = int((kmax-kmin)/dk)

	integrate_epsilon  = 0.d0
	umag0 = 2 * vis * kmin**2 * energy_init(kmin,epsf,l_e,eta,c_eta,c_l)
	do i=1,ikrange
		kmode             = kmin + i * dk
		umag              = 2 * vis * kmode**2 * energy_init(kmode,epsf,l_e,eta,c_eta,c_l)
		integrate_epsilon = integrate_epsilon + (umag+umag0) * dk / 2
		umag0             = umag
    end do
end function integrate_epsilon




!^^^^^^^^^^^ Subroutines related to the determination of C_eta and C_l
subroutine xmnewt(c_eta,c_l,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	real(prcn),intent(inout) :: c_eta,c_l
	real(prcn),intent(in)    :: tke,eps,l_e,eta,kmax,vis
	INTEGER NTRIAL,N,NP
	real(prcn) TOLX,TOLF
	PARAMETER(NTRIAL=50,TOLX=1.0D-10,N=2,TOLF=1.0D-10,NP=2)
	INTEGER i,j,k,kk
	real(prcn) xx,fjac(NP,NP),fvec(NP),x(NP)

	x(1)=c_eta
	x(2)=c_l

	call usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
!	write(*,'(/1x,t5,a,t14,a,t29,a/)') 'I','X(I)','F'
!	do i=1,N
!		write(*,'(1x,i4,2e15.6)') i,x(i),fvec(i)
!	enddo

	do j=1,NTRIAL
		call mnewt(1,x,N,TOLX,TOLF,tke,eps,l_e,eta,kmax,vis)
		call usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
		write(*,'(/1x,t5,a,t14,a,t29,a/)') 'I','X(I)','F'
		do i=1,N
			write(*,'(1x,i4,2e15.6)') i,x(i),fvec(i)
		enddo
		if (abs(fvec(1))<1D-10.and.abs(fvec(2))<1D-10) exit
	enddo
	c_eta = x(1)
	c_l   = x(2)
end subroutine xmnewt


SUBROUTINE usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER i,n
	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) fjac(n,n),fvec(n),x(n)

	call funcv(n,x,fvec,tke,eps,l_e,eta,kmax,vis)
	call fdjac(n,x,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
END

SUBROUTINE funcv(n,x,f,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER n
	real(prcn) x(n),f(n)

	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) :: c_eta,c_l
	real(prcn) :: krange,kmode,dk,sum1,sum2,umag0,umag
	integer :: i,ikrange
	real(prcn), external  :: energy_init,integrate_energy,integrate_epsilon

	c_eta = x(1)
	c_l   = x(2)
	
	!tke: integral of E(k): 0-kmax
	sum1=integrate_energy(1d-6,kmax,eps,l_e,eta,c_eta,c_l)
	f(1)=tke-sum1

	!eps: integratal of 2*nu*k^2*E(k)*dk (0-kmax)
	sum2=integrate_epsilon(1d-6,kmax,vis,eps,l_e,eta,c_eta,c_l)
	f(2)=eps-sum2
	return
END

SUBROUTINE fdjac(n,x,fvec,df,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER n
	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) df(n,n),fvec(n),x(n)
	real(prcn) small
	PARAMETER (small=1.e-4)
	INTEGER i,j
	real(prcn) h,temp,f(n)

	do j=1,n
		temp=x(j)
		h=small*abs(temp)
		if(h.eq.0.)h=small
		x(j)=temp+h
		h=x(j)-temp
		call funcv(n,x,f,tke,eps,l_e,eta,kmax,vis)
		x(j)=temp
		do i=1,n
			df(i,j)=(f(i)-fvec(i))/h
		enddo
	enddo
	return
END

SUBROUTINE mnewt(ntrial,x,n,tolx,tolf,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER n,ntrial
	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) tolf,tolx,x(n)
	INTEGER i,k,indx(n)
	real(prcn) d,errf,errx,fjac(n,n),fvec(n),p(n)

	do k=1,ntrial
		call usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
		errf=0.
		do i=1,n
			errf=errf+abs(fvec(i))
		enddo
		if(errf.le.tolf)return
		do i=1,n
			p(i)=-fvec(i)
		enddo
		call ludcmp(fjac,n,indx,d)
		call lubksb(fjac,n,indx,p)
		errx=0.
		do i=1,n
			errx=errx+abs(p(i))
			if (x(i)+p(i)>0) then
				x(i)=x(i)+p(i)
			else
				x(i)=x(i)*0.5
			endif
		enddo
		if(errx.le.tolx)return
	enddo
	return
END


SUBROUTINE lubksb(a,n,indx,b)
	use precision
	implicit none
	INTEGER n
	integer indx(n)
	real(prcn) a(n,n),b(n)
	INTEGER i,ii,j,ll
	real(prcn) sum
	ii=0
	do i=1,n
		ll=indx(i)
		sum=b(ll)
		b(ll)=b(i)
		if (ii.ne.0)then
			do j=ii,i-1
				sum=sum-a(i,j)*b(j)
			enddo
		else if (sum.ne.0.) then
			ii=i
		endif
		b(i)=sum
	enddo
	do i=n,1,-1
		sum=b(i)
		do j=i+1,n
			sum=sum-a(i,j)*b(j)
		enddo
		b(i)=sum/a(i,i)
	enddo
	return
END

SUBROUTINE ludcmp(a,n,indx,d)
	use precision
	implicit none
	INTEGER n
	integer indx(n),NMAX
	real(prcn) d,a(n,n),TINY
	PARAMETER (NMAX=500,TINY=1.0e-20)
	INTEGER i,imax,j,k
	real(prcn) aamax,dum,sum,vv(n)
	d=1.
	do i=1,n
		aamax=0.
		do j=1,n
			if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
		enddo
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
	enddo
	do j=1,n
		do i=1,j-1
			sum=a(i,j)
			do k=1,i-1
				sum=sum-a(i,k)*a(k,j)
			enddo
			a(i,j)=sum
		enddo
		aamax=0.
		do i=j,n
			sum=a(i,j)
			do k=1,j-1
				sum=sum-a(i,k)*a(k,j)
			enddo
			a(i,j)=sum
			dum=vv(i)*abs(sum)
			if (dum.ge.aamax) then
				imax=i
				aamax=dum
			endif
		enddo
		if (j.ne.imax)then
			do k=1,n
				dum=a(imax,k)
				a(imax,k)=a(j,k)
				a(j,k)=dum
			enddo
			d=-d
			vv(imax)=vv(j)
		endif
		indx(j)=imax
		if(a(j,j).eq.0.)a(j,j)=TINY
		if(j.ne.n)then
			dum=1./a(j,j)
			do i=j+1,n
				a(i,j)=a(i,j)*dum
			enddo
		endif
	enddo
	return
END
!--------------------------------------------------

!^^^^^^^^^^^ Subroutines related to the mean value theorem
subroutine xrtnewt(root,k0,kmax,epsf,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) epsf,l_e,eta,c_eta,c_l
	INTEGER N,NBMAX
	real(prcn):: X1,X2
	PARAMETER(N=100,NBMAX=1)
	INTEGER i,nb
	real(prcn):: integrate_energy,f,rtnewt,root,xacc,xb1(NBMAX),xb2(NBMAX)
	EXTERNAL funcd,f,integrate_energy

	real(prcn) :: dk,sum,umag0,umag,kmode,k0,kmax,const
	integer :: ikrange

	const = integrate_energy(k0,kmax,epsf,l_e,eta,c_eta,c_l)
	const = const / (kmax-k0)

	X1=k0
	X2=kmax
	nb=NBMAX
	call zbrak(f,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,X1,X2,N,xb1,xb2,nb)
!	write(*,'(/1x,a)') 'Roots of BESSJ0:'
!	write(*,'(/1x,t19,a,t31,a/)') 'x','F(x)'
	do i=1,nb
		xacc=(1.0e-6)
		root=rtnewt(funcd,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,xb1(i),xb2(i),xacc)
!		write(*,'(1x,a,i2,2x,f12.6,e16.4)') 'Root ',i,root,f(root,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	enddo
END

function f(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) :: x,f
	real(prcn) :: k0,kmax,const,tke,epsf,l_e,eta,c_eta,c_l
	real(prcn),external :: energy_init

	f=energy_init(x,epsf,l_e,eta,c_eta,c_l)-const
end

function df(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,f)
	use precision
	implicit none
	real(prcn)  :: x,f,df
	real(prcn)  :: k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn),external  :: fdjac1
	
!	df=2*x
	df=fdjac1(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,f)
end

function fdjac1(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,fvec)
	use precision
	implicit none
	INTEGER n
	real(prcn) fdjac1
	real(prcn) x,fvec,df
	real(prcn) k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn) small
	PARAMETER (small=1.e-4)
	INTEGER i,j
	real(prcn) h,temp
	real(prcn),external :: f

	temp=x
	h=small*abs(temp)
	if(h.eq.0.)h=small
	x=temp+h
	h=x-temp

	fdjac1=(f(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)-fvec)/h
	x=temp
	return
END

SUBROUTINE funcd(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,fn,dfn)
	use precision
	implicit none
	real(prcn)  :: fn,dfn,x
	real(prcn)  :: k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn),external :: f,df

	fn  =  f(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	dfn = df(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,fn)
	return
END

FUNCTION rtnewt(funcd,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,x1,x2,xacc)
	use precision
	implicit none
	INTEGER JMAX
	real(prcn):: rtnewt,x1,x2,xacc
	real(prcn)::k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	EXTERNAL funcd
	PARAMETER (JMAX=200)
	INTEGER j
	real(prcn):: df,dx,f
	rtnewt=.5*(x1+x2)
	do j=1,JMAX
		call funcd(rtnewt,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,f,df)
		dx=f/df
		rtnewt=rtnewt-dx
!		if((x1-rtnewt)*(rtnewt-x2).lt.0.) pause
!		if((x1-rtnewt)*(rtnewt-x2).lt.0.) write (*,*) 'rtnewt jumped out of brackets'
		if(abs(dx).lt.xacc) return
	enddo
!	write (*,*) 'rtnewt exceeded maximum iterations'
!	pause 
END

SUBROUTINE zbrak(fx,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,x1,x2,n,xb1,xb2,nb)
	use precision
	implicit none
	INTEGER n,nb
	real(prcn):: x1,x2,xb1(nb),xb2(nb)
	real(prcn) :: k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn),EXTERNAL :: fx
	INTEGER i,nbb
	real(prcn):: dx,fc,fp,x
	nbb=0
	x=x1
	dx=(x2-x1)/n
	fp=fx(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	do i=1,n
		x=x+dx
		fc=fx(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
		if(fc*fp.le.0.) then
			nbb=nbb+1
			xb1(nbb)=x-dx
			xb2(nbb)=x
			if(nbb.eq.nb) exit
		endif
		fp=fc
	enddo
	nb=nbb
	return
END


#if 0
		!for i=1,j=1,k<>1 line
		if (local_o_start(1)==1 .and. local_o_start(2)==1) then
			i=1
			j=1
			do k=global_n(3)/2+1,global_n(3)
				u(k,i,j,1) = conjg(u(global_n(3)+2-k,i,j,1))
				u(k,i,j,2) = conjg(u(global_n(3)+2-k,i,j,2))
				u(k,i,j,3) = conjg(u(global_n(3)+2-k,i,j,3))
			enddo
		endif

		!for i=1, j<>1, k=1 line
		if (local_o_start(1)==1 .and. local_o_start(3)==1) then
			if (nprocy==1) then
				i=1
				k=1
				do j=2,global_n(2)/2
					u(k,i,global_n(2)+2-j,1) = dcmplx(dble(u(k,i,j,1)),(-1.0)*aimag(u(k,i,j,1)))
					u(k,i,global_n(2)+2-j,2) = dcmplx(dble(u(k,i,j,2)),(-1.0)*aimag(u(k,i,j,2)))
					u(k,i,global_n(2)+2-j,3) = dcmplx(dble(u(k,i,j,3)),(-1.0)*aimag(u(k,i,j,3)))
				enddo
			else
			endif
		endif

		!for i=1, J<>1, k<>1 
		if (local_o_start(1)==1) then
			if (nprocy==1) then
				i=1
				do k=2, global_n(3)/2
					do j=2, global_n(2)/2
						u(global_n(3)+2-k,i,global_n(2)+2-j,:) = conjg(u(k,              i,j,:))
						u(k,              i,global_n(2)+2-j,:) = conjg(u(global_n(3)+2-k,i,j,:))
					end do
				enddo
			else

			endif
		endif
#endif		




