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

module physalis_mod
#include "../FLO/ibm.h"
	use global_data
	use nlmainarrays , only : ubcp, pbcp
!	use fftw_interface
	use constants
	use nlarrays
!	use randomno
	use general_funcs
!	use string_funcs
!	use init_turb , only : avr, tke, iauto_on, autocorrelation
	use boundary_condition
	use dependent_functions
	use bcsetarrays
!	use epsarray
	use postproc_funcs
	use nr
	use nrtype
	use nrutil
	use mypost_process
	implicit none

contains

subroutine physalis
	implicit none

	type :: node_type
		integer :: ind
		real(8) :: pos(ndim)
	endtype node_type

	type :: cage_type
		integer :: num
		type(node_type), allocatable :: node(:)
	endtype cage_type

	type(cage_type), allocatable, dimension(:) :: cage_pres, cage_vel

	integer :: i, j, k, l, i1, j1, k1, ii, jj, kk, iii, jjj, kkk, m, n, o, ijk, ijk_1
	integer :: im, ip, jm, jp, km, kp, ip1, jp1, kp1, i2, j2, k2, idim
	
	real(8) :: tmp_vec1(ndim), tmp_vec2(ndim), cart(ndim), sphr(ndim), dist1, dist2, dist3, dist4
	integer :: tmp_num, exclude_num
	integer, allocatable :: tmp_ind(:), exclude_ind(:)
	logical :: neighb_insolid, registered, select_node


	real(8), allocatable, dimension(:,:) :: p_nm, pt_nm, f_nm,ft_nm, x_nm, xt_nm, coeff, u_mat, v_mat
	real(8), allocatable, dimension(:)   :: coeff_tmp, rhs, w_mat, x_mat, rhs_check
	integer :: nc, np, nv, iv, ibody, mm, nn, eq_num
	real(8) :: rad, theta, phi, s
	real(8) :: tmp1, tmp2, tmp3, w_min, w_max, error, denom
	character*50 filename

!---------------------------------------------------
	nc = 4

	call form_cage_nodes
!	call cage_output
	call change_coordinates
	call exact_field

	contains

	subroutine exact_field
		implicit none

		integer :: iphs
		real(8) :: da(2), mean_pres, pres, p_th, vel_th(ndim), dpdr_th, dp_rdtheta_th, dp_rsindphi_th, d2urdr2_th ,d2uthetadr2_th, d2uphidr2_th, dissip_pres, dissip_visc
		real(8), dimension(ndim):: vort, xl, xpb, pres_total1, visc_total1, is, ul, ppll
		integer, dimension(ndim):: pcellb, vcellb
		real(8) :: pl, pres_drag1, visc_drag1, total_drag1

		write (*,*) "SOLVING FOR THE EXACT FIELD"


		allocate(p_nm(nc,0:nc),pt_nm(nc,0:nc))
		allocate(f_nm(nc,0:nc),ft_nm(nc,0:nc))
		allocate(x_nm(nc,0:nc),xt_nm(nc,0:nc))

		! OBTAINING THE COEFFICIENTS Pnm, P~nm, PHInm, PHI~nm, Xnm, X~nm
		do ibody=1, nbody
			np = cage_pres(ibody)%num
!			np = 0
			nv = cage_vel(ibody)%num
			mm = np + nv*3
!			mm = nv
			nn = 6 * (nc)*(nc+1)
!			allocate(p_nm(nc,0:nc),pt_nm(nc,0:nc))
!			allocate(f_nm(nc,0:nc),ft_nm(nc,0:nc))
!			allocate(x_nm(nc,0:nc),xt_nm(nc,0:nc))

			p_nm=0d0; pt_nm=0d0
			f_nm=0d0; ft_nm=0d0
			x_nm=0d0; xt_nm=0d0

			allocate(coeff(mm,nn), rhs(mm))
			coeff = 0d0
			rhs   = 0d0

			allocate(coeff_tmp(nn))
			coeff_tmp = 0d0

			eq_num = 0 ! counter for equations
			! THE EQUATION FOR EACH PRESSURE CAGE NODE

			do ip=1, np
				eq_num = eq_num + 1

				sphr(:) = cage_pres(ibody)%node(ip)%pos(:)
				rad   = sphr(1)
				theta = sphr(2)
				phi   = sphr(3)
				s     = rad/radbdy(ibody)
!				s     = radbdy(ibody)*dy

				l=0 ! counter for coefficients
				do n=1, nc
					tmp1 = (s**n - n*(2d0*n-1d0)/(2d0*(n+1d0))*s**(-n-1d0)) * (rhof*vis**2/(radbdy(ibody)*dy)**2)
					tmp2 = (-n*(4d0*n**2-1d0)/(n+1d0)*s**(-n-1d0))          * (rhof*vis**2/(radbdy(ibody)*dy)**2)
					! Pnm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * dcos(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! P~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * dsin(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! PHInm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * dcos(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! PHI~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * dsin(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! Xnm
					do m=0, nc
						l = l+1
						coeff_tmp(l) = 0d0
					enddo
					! X~nm
					do m=0, nc
						l = l+1
						coeff_tmp(l) = 0d0
					enddo
				enddo
				call sphr_to_cart(sphr,ibody,cart)
				cart(:) = cart(:)*dy!-doml(:)/2
				mean_pres = dot_product(mpg,cart)

				ijk = cage_pres(ibody)%node(ip)%ind
				call ind1t3(ijk,i,j,k)
				coeff(eq_num,:) = coeff_tmp(:)
				rhs(eq_num) = pbcp(i,j,k)+mean_pres  !P_0 is not added to the RHS, ?

				! CHECK POINT
				if (l/=nn) then
					write (*,*) "THE NUMBER OF COEFFICIENTS DOES NOT MATCH THE PREDICTION"
					write (*,*) "STOP"
					stop
				endif
			enddo

			! THE EQUATION FOR VELOCITY CAGE NODE U_r
			do iv=1, nv
				eq_num = eq_num + 1

				sphr(:) = cage_vel(ibody)%node(iv)%pos(:)
				rad   = sphr(1)
				theta = sphr(2)
				phi   = sphr(3)
				s = rad/radbdy(ibody)

				l=0 ! counter for coefficients
				do n=1, nc
					tmp1 = (n/2d0/(2d0*n+3d0)*s**(n+1d0) - n/4d0*s**(-n) + n*(2d0*n+1d0)/4d0/(2d0*n+3d0)*s**(-n-2d0)) * vis/(radbdy(ibody)*dy)
					tmp2 = (n*s**(n-1d0) - n*(2d0*n+1d0)/2d0*s**(-n) + n*(2d0*n-1d0)/2d0*s**(-n-2d0))                 * vis/(radbdy(ibody)*dy)
					! Pnm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * dcos(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! P~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * dsin(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! PHInm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * dcos(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! PHI~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * dsin(m*phi) * plgndr(n,m,dcos(theta))
						endif
					enddo
					! Xnm
					do m=0, nc
						l = l+1
						coeff_tmp(l) = 0d0
					enddo
					! X~nm
					do m=0, nc
						l = l+1
						coeff_tmp(l) = 0d0
					enddo
				enddo
				ijk = cage_vel(ibody)%node(iv)%ind
				call ind1t3(ijk,i,j,k)
				sphr(:) = ubcp(i,j,k,:) !+ufmean(:)

				call vec_cart_to_sphr(sphr,theta,phi)

				coeff(eq_num,:) = coeff_tmp(:)
				rhs(eq_num) = sphr(1)

				! CHECK POINT
				if (l/=nn) then
					write (*,*) "THE NUMBER OF COEFFICIENTS DOES NOT MATCH THE PREDICTION"
					write (*,*) "STOP"
					stop
				endif
			enddo


			! THE EQUATION FOR VELOCITY CAGE NODE U_theta
			do iv=1, nv
				eq_num = eq_num + 1

				sphr(:) = cage_vel(ibody)%node(iv)%pos(:)
				rad   = sphr(1)
				theta = sphr(2)
				phi   = sphr(3)
				s = rad/radbdy(ibody)
				l=0 ! counter for coefficients
				do n=1, nc
					tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)*s**(n+1d0)+(n-2d0)/4d0/(n+1d0)*s**(-n)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)*s**(-n-2d0)) &
						& * vis/(radbdy(ibody)*dy)
					tmp2 = (s**(n-1d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*s**(-n)-n*(2d0*n-1d0)/2d0/(n+1d0)*s**(-n-2d0)) * vis/(radbdy(ibody)*dy)
					tmp3 = (s**n-s**(-n-1d0)) * vis/(radbdy(ibody)*dy)
					! Pnm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * dcos(m*phi) * dplgndr(n,m,theta)
						endif
					enddo

					! P~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * dsin(m*phi) * dplgndr(n,m,theta)
						endif
					enddo
					! PHInm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * dcos(m*phi) * dplgndr(n,m,theta)
						endif
					enddo
					! PHI~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * dsin(m*phi) * dplgndr(n,m,theta)
						endif
					enddo
					! Xnm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp3 * (-m) * dsin(m*phi) * plgndr(n,m,dcos(theta)) / dsin(theta)
						endif
					enddo
					! X~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp3 * m * dcos(m*phi) * plgndr(n,m,dcos(theta)) / dsin(theta)
						endif
					enddo
				enddo
				ijk = cage_vel(ibody)%node(iv)%ind
				call ind1t3(ijk,i,j,k)
				sphr(:) = ubcp(i,j,k,:) !+ufmean(:)
				call vec_cart_to_sphr(sphr,theta,phi)

				coeff(eq_num,:) = coeff_tmp(:)
				rhs(eq_num) = sphr(2)

				! CHECK POINT
				if (l/=nn) then
					write (*,*) "THE NUMBER OF COEFFICIENTS DOES NOT MATCH THE PREDICTION"
					write (*,*) "STOP"
					stop
				endif
			enddo


			! THE EQUATION FOR VELOCITY CAGE NODE U_phi
			do iv=1, nv
				eq_num = eq_num + 1

				sphr(:) = cage_vel(ibody)%node(iv)%pos(:)
				rad   = sphr(1)
				theta = sphr(2)
				phi   = sphr(3)
				s = rad/radbdy(ibody)
				l=0 ! counter for coefficients
				do n=1, nc
					tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)*s**(n+1d0)+(n-2d0)/4d0/(n+1d0)*s**(-n)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)*s**(-n-2d0)) &
						& * vis/(radbdy(ibody)*dy)
					tmp2 = (s**(n-1d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*s**(-n)-n*(2d0*n-1d0)/2d0/(n+1d0)*s**(-n-2d0)) * vis/(radbdy(ibody)*dy)
					tmp3 =-(s**n-s**(-n-1d0)) * vis/(radbdy(ibody)*dy)
					! Pnm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * (-m) * dsin(m*phi) * plgndr(n,m,dcos(theta)) / dsin(theta)
						endif
					enddo

					! P~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp1 * m * dcos(m*phi) * plgndr(n,m,dcos(theta)) / dsin(theta)
						endif
					enddo

					! PHInm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * (-m) * dsin(m*phi) * plgndr(n,m,dcos(theta)) / dsin(theta)
						endif
					enddo

					! PHI~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp2 * m * dcos(m*phi) * plgndr(n,m,dcos(theta)) / dsin(theta)
						endif
					enddo

					! Xnm
					do m=0, nc
						l = l+1
						if (m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp3 * dcos(m*phi) * dplgndr(n,m,theta)
						endif
					enddo

					! X~nm
					do m=0, nc
						l = l+1
						if (m==0.or.m>n) then
							coeff_tmp(l) = 0d0
						else
							coeff_tmp(l) = tmp3 * dsin(m*phi) * dplgndr(n,m,theta)
						endif
					enddo
				enddo
				ijk = cage_vel(ibody)%node(iv)%ind
				call ind1t3(ijk,i,j,k)
				sphr(:) = ubcp(i,j,k,:) !+ufmean(:)
				call vec_cart_to_sphr(sphr,theta,phi)

				coeff(eq_num,:) = coeff_tmp(:)
				rhs(eq_num) = sphr(3)

				! CHECK POINT
				if (l/=nn) then
					write (*,*) "THE NUMBER OF COEFFICIENTS DOES NOT MATCH THE PREDICTION"
					write (*,*) "STOP"
					stop
				endif
			enddo

			! CHECK POINT
			if (eq_num/=mm) then
				write (*,*) "THE NUMBER OF EQUATIONS DOES NOT MATCH THE PREDICTION"
				write (*,*) "STOP"
				stop
			endif

			write (*,*) "CALLING THE LIST SQUARE ROUTINE"
!			allocate(u_mat(mm,nn),w_mat(nn),v_mat(nn,nn),x_mat(nn),rhs_check(mm))
			allocate(u_mat(mm,nn),x_mat(nn),rhs_check(mm))
			u_mat(1:mm,1:nn)=coeff(1:mm,1:nn)
			! decompose matrix u_mat
!			call svdcmp(u_mat(1:mm,1:nn),w_mat(1:nn),v_mat(1:nn,1:nn))
			! find maximum singular value	subroutine p_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
!			w_max=max(maxval(w_mat(1:nn)),0d0)
			! define "small"
!			w_min=w_max*(1.0d-10)
			! zero the "small" singular values
!			where (w_mat(1:nn) < w_min) w_mat(1:nn)=0d0
			! backsubstitute for each right-hand side vector
!			call svbksb(u_mat(1:mm,1:nn),w_mat(1:nn),v_mat(1:nn,1:nn),rhs(1:mm),x_mat(1:nn))
			call xwnnls(u_mat(1:mm,1:nn),rhs(1:mm),x_mat(1:nn),mm,nn)
			rhs_check(1:mm)=matmul(coeff(1:mm,1:nn),x_mat(1:nn))
			error = 0.0
			do i=1, mm
				error = error + abs(rhs_check(i)-rhs(i))**2
			enddo
			denom = 0.0
			do i=1, mm
				denom = denom + abs(rhs(i))**2
			enddo
			error = sqrt(error/denom) !error/denom
			write(*,'(1A,1D15.7)') "THE ERROR BETWEEN ACTUAL AND ESTIMATED RHS = ", error

!			open(unit=1, file="RHS_N4_post.dat", status="replace", action="write")
!			write (1,"(2F15.7)") ((rhs(i),rhs_check(i)), i=1,mm)
!			close (1)
!			open(unit=1, file="VEC_N4_post.dat", status="replace", action="write")
!			write (1,"(1F15.7)") ((x_mat(i)), i=1,nn)
!			close (1)

			l=0
			do n=1, nc
				! Pnm
				do m=0, nc
					l = l+1
					p_nm(n,m) = x_mat(l)
				enddo

				! P~nm
				do m=0, nc
					l = l+1
					pt_nm(n,m) = x_mat(l)
				enddo

				! PHInm
				do m=0, nc
					l = l+1
					f_nm(n,m) = x_mat(l)
				enddo

				! PHI~nm
				do m=0, nc
					l = l+1
					ft_nm(n,m) = x_mat(l)
				enddo

				! Xnm
				do m=0, nc
					l = l+1
					x_nm(n,m) = x_mat(l)
				enddo

				! X~nm
				do m=0, nc
					l = l+1
					xt_nm(n,m) = x_mat(l)
				enddo
			enddo

			! CHECK POINT
			if (l/=nn) then
				write (*,*) "THE NUMBER OF COEFFICIENTS DOES NOT MATCH THE PREDICTION"
				write (*,*) "STOP"
				stop
			endif

			open(unit=1, file="P_nm_N4_post.dat", status="replace", action="write")
			write (1,"(5F15.7)") ((p_nm(i,j), j=0,nc), i=1,nc)
			close (1)

			open(unit=1, file="Pt_nm_N4_post.dat", status="replace", action="write")
			write (1,"(5F15.7)") ((pt_nm(i,j), j=0,nc), i=1,nc)
			close (1)
!
			open(unit=1, file="F_nm_N4_post.dat", status="replace", action="write")
			write (1,"(5F15.7)") ((f_nm(i,j), j=0,nc), i=1,nc)
			close (1)

			open(unit=1, file="Ft_nm_N4_post.dat", status="replace", action="write")
			write (1,"(5F15.7)") ((Ft_nm(i,j), j=0,nc), i=1,nc)
			close (1)

			open(unit=1, file="X_nm_N4_post.dat", status="replace", action="write")
			write (1,"(5F15.7)") ((x_nm(i,j), j=0,nc), i=1,nc)
			close (1)

			open(unit=1, file="Xt_nm_N4_post.dat", status="replace", action="write")
			write (1,"(5F15.7)") ((xt_nm(i,j), j=0,nc), i=1,nc)
			close (1)
!------------------------------------------------
			iphs = 1!part_array(ibody)%iphs
			nbnd = phase_array(iphs)%nbnd
			NULLIFY(bndarray)
			bndarray => phase_array(iphs)%bndpts
			da(1)=4.*pi*(radbdy(ibody)*dx)**2./real(nbnd,prcn)

			open(unit=1, file="vorticity_exect_N4_post.dat")

			pres_total1 = 0d0
			visc_total1 = 0d0
			do l=1,nbnd
!				rad = zero
!				do n=1,ndim
!					xl(n)=xc(ibody,n)+ bndarray(n,l)*radbdy(ibody)
!					rad=rad+(bndarray(n,l)*radbdy(ibody))**2.0
!				enddo
!				rad = dsqrt(rad)


!				xpb(1) = xl(1)-0.5
!				xpb(2:3)=xl(2:3)
!				do n = 1, ndim
!					if(xpb(n).lt.zero) then 
!						pcellb(n) = int(xpb(n)-1)
!					else 
!						pcellb(n) = int(xpb(n))
!					endif
!					if(xl(n).lt.zero) then 
!						vcellb(n) = int(xl(n)-1)
!					else 
!						vcellb(n) = int(xl(n))
!					endif
!				enddo
!#if PARALLEL
!				xltemp = xl(1)
!				xptemp = xpb(1)
!				pcelltemp = pcellb(1)
!				vcelltemp = vcellb(1)
!				if(l.eq.FOCUS_POINT)then
!					PRINT*,' xl = ', myid, xltemp, xptemp, pcelltemp
!				endif
!				if(.not.CELL_IN_VEL_GRID(vcelltemp))then
!					WEST_PERIODIC_IMAGE(vcellb(1),vcelltemp,xl(1),xltemp)
!					WEST_PERIODIC_IMAGE(pcellb(1),pcelltemp,xpb(1),xptemp)
!					EAST_PERIODIC_IMAGE( vcellb(1),vcelltemp,xl(1),xltemp)
!					EAST_PERIODIC_IMAGE_PRES(pcellb(1),pcelltemp,xpb(1),xptemp)
!					if(l.eq.FOCUS_POINT)then
!						PRINT*,' xl IMAGES = ', myid,xltemp, xptemp, pcelltemp
!					endif
!
!!					if(.not.CELL_IN_VEL_GRID(vcelltemp)) goto 777
!				endif
!				vcellb(1) = vcelltemp
!				pcellb(1) = pcelltemp
!				xl(1) = xltemp
!				xpb(1) = xptemp
!#endif
!				pl=zero
!				call interpolate_pdata(pcellb,xpb,ppll,pl,l)
!				call interpolate_udata(vcellb,xl,ib&
!				  &,ie,jb,je,kb,ke,ul,nll,onll,dfll, 0,m, l, onew) 
		  
				cart(:) = bndarray(:,l)*radbdy(ibody)
				call cart_to_sphr(cart,ibody,sphr)
				call p_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
				pres_total1(:) = pres_total1(:) + pres*cd(:,l)*da(1)

				if (abs(sphr(2))>1d-3.and.abs(sphr(2)-pi)>1d-3) then
					call omega_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,vort)
!					call vel_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,vort)
					tmp_vec1(1) = 0d0
					tmp_vec1(2) = vort(3)
					tmp_vec1(3) =-vort(2)
					call vec_sphr_to_cat(tmp_vec1,sphr(2),sphr(3))
!					tmp_vec1 = -tmp_vec1
					visc_total1(:) = visc_total1(:) + tmp_vec1(:)*vis*da(1)		

					write (1,"(6D15.7)") cart(:),tmp_vec1(:)
				else
					write (1,"(6D15.7)") cart(:),0d0,0d0,0d0
				endif
	 		enddo
			close (1)

			write (*,"(A,3D15.7)") "VISC_PART =",visc_total1(:)
			write (*,"(A,3D15.7)") "pres_PART =",pres_total1(:)


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			visc_drag1 = DSQRT(dot_product(visc_total1(1:ndim), visc_total1(1:ndim)))
			pres_drag1 = DSQRT(dot_product(pres_total1(1:ndim), pres_total1(1:ndim)))

			pres_drag1 = pres_drag1/norm_factor/(one-maxvolfrac)
			visc_drag1 = visc_drag1/norm_factor/(one-maxvolfrac)
			total_drag1 = (pres_drag1+visc_drag1)

!			open (unit=1,file="drag_comp_N4_post.dat",status="pld",action="write",position="append")
			open (unit=1,file="drag_comp_N4_post.dat",status="replace",action="write")
			write (1,"(A)") "variables=D<sub>m</sub>,<greek>f</greek>,Visc.,Pres.,Total"
			write (1,"(A)") "zone"
			write (1,"(1I8,1F8.2,3D15.7)") int(dbydx),maxvolfrac,visc_drag,pres_drag,total_drag
			write (1,"(A)") "zone"
			write (1,"(1I8,1F8.2,3D15.7)") int(dbydx),maxvolfrac,visc_drag1,pres_drag1,total_drag1
			close (1)
!-------------------------------


!tmp_vec1(1) = pi*vis**2*(p_nm(1,1)+6*f_nm(1,1))
!tmp_vec1(2) = pi*vis**2*(pt_nm(1,1)+6*ft_nm(1,1))
!tmp_vec1(3) = pi*vis**2*(p_nm(1,0)+6*f_nm(1,0))

!write (*,"(3D15.7)") tmp_vec1(:)
!write (*,"(3D15.7)") tmp_vec1(:)/(3*pi*radbdy(ibody)*dy*vis*umeanslip*(1.-maxvolfrac))

!			call reynolds_stress_tensor
!			tke = tke * half*dot_product(ufmean(:)-usmean(:),ufmean(:)-usmean(:))

			filename = trim(run_name)//"_dissip_physalis_N4.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")		

			write (1,"(A)") "variables=<greek>q</greek>,p/<greek>r</greek>|<math>a</math>W<math>q</math>|<sup>2</sup>,&
						  &u<sub>r0</sub>,u<sub><greek>q</greek>0</sub>,u<sub><greek>f</greek>0</sub>,&
						  &u<sub>r1</sub>,u<sub><greek>q</greek>1</sub>,u<sub><greek>f</greek>1</sub>,&
						  &<math>6</math>p/<math>6</math>r,<math>6</math>p/r<math>6</math><greek>q</greek>,&
						  &<math>6</math>p/rsin(<greek>q</greek>)<math>6</math><greek>f</greek>,&
						  &(<math>Q</math><sup>2</sup><b>U</b>)<sub>r</sub>,&
						  &(<math>Q</math><sup>2</sup><b>U</b>)<sub><greek>q</greek></sub>,&
						  &(<math>Q</math><sup>2</sup><b>U</b>)<sub><greek>f</greek></sub>,&
						  &<greek>Q</greek><sub>pres</sub>,<greek>Q</greek><sub>visc</sub>,&
						  &<greek>Q</greek><sub>total</sub>"

			do i=1, 359
				p_th    = 0d0
				vel_th  = 0d0
				dpdr_th        = 0d0
				dp_rdtheta_th  = 0d0
				dp_rsindphi_th = 0d0
				d2urdr2_th     = 0d0
				d2uthetadr2_th = 0d0
				d2uphidr2_th   = 0d0

				do j=1, 360
					rad = radbdy(ibody)
					sphr(1) = rad
					sphr(2) = i*pi/180.0/2.0
					sphr(3) = j*pi/180.0
					call p_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
					p_th = p_th + pres
					call dpdr_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
					dpdr_th = dpdr_th + pres
					call dp_rdtheta_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
					dp_rdtheta_th = dp_rdtheta_th + pres
					call dp_rsindphi_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
					dp_rsindphi_th = dp_rsindphi_th + pres

					call d2urdr2_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,tmp1)
					d2urdr2_th = d2urdr2_th + tmp1
					call d2uthetadr2_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,tmp1)
					d2uthetadr2_th = d2uthetadr2_th + tmp1
					call d2uphidr2_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,tmp1)
					d2uphidr2_th = d2uphidr2_th + tmp1


					sphr(1) = rad + dr
					call vel_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,tmp_vec1)
					vel_th(:) = vel_th(:) + tmp_vec1(:)
				enddo
				tmp_vec2 = -(ufmean-usmean)
				call vec_cart_to_sphr(tmp_vec2,sphr(2),0d0)

				p_th   = p_th   /360
				vel_th = vel_th /360

				dpdr_th        =  dpdr_th        /360
				dp_rdtheta_th  =  (-1.) * dp_rdtheta_th  /360       ! probably an error (a minus sign) here
				dp_rsindphi_th =  dp_rsindphi_th /360

				d2urdr2_th     = vis*d2urdr2_th     /360
				d2uthetadr2_th = vis*d2uthetadr2_th /360
				d2uphidr2_th   = vis*d2uphidr2_th   /360

				tmp_vec1(1) = dpdr_th
				tmp_vec1(2) = dp_rdtheta_th
				tmp_vec1(3) = dp_rsindphi_th
				dissip_pres = dot_product(tmp_vec1,tmp_vec2)

				tmp_vec1(1) = d2urdr2_th
				tmp_vec1(2) = d2uthetadr2_th
				tmp_vec1(3) = d2urdr2_th
				dissip_visc = dot_product(tmp_vec1,tmp_vec2)

				p_th     = p_th     / (rhof*umeanslip**2)
				vel_th   = vel_th   / umeanslip
				tmp_vec2 = tmp_vec2 / umeanslip

				dpdr_th        = dpdr_th        / (rhof*umeanslip**2/dia_phys)
				dp_rdtheta_th  = dp_rdtheta_th  / (rhof*umeanslip**2/dia_phys)
				dp_rsindphi_th = dp_rsindphi_th / (rhof*umeanslip**2/dia_phys)

				d2urdr2_th     = d2urdr2_th     / (umeanslip/dia_phys**2)
				d2uthetadr2_th = d2uthetadr2_th / (umeanslip/dia_phys**2)
				d2uphidr2_th   = d2uphidr2_th   / (umeanslip/dia_phys**2)

				dissip_pres = dissip_pres / (vis*(umeanslip/dia_phys)**2)
				dissip_visc = dissip_visc / (vis*(umeanslip/dia_phys)**2)
				write (1,"(17D15.7)") sphr(2)*180.0/pi, p_th, tmp_vec2(:), vel_th(:), dpdr_th, dp_rdtheta_th, dp_rsindphi_th,&
								& d2urdr2_th, d2uthetadr2_th, d2uphidr2_th, dissip_pres, dissip_visc, dissip_pres+dissip_visc
			enddo
write (*,*) umeanslip
			close (1)
#if 0
			do j=0, 180
				write (1,*) "zone"
				do i=1, 359
					rad = radbdy(ibody)
					sphr(1) = rad
					sphr(2) = i*pi/180.0/2.0
					sphr(3) = j*pi/180.0*2
					call p_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)

					sphr(1) = rad + dr
					call vel_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,tmp_vec1)
					write (1,"(5D15.7)") sphr(2)*180.0/pi, pres, tmp_vec1(:)
				enddo
			enddo
			close (1)

			do iv=1,nv
				sphr(:) = cage_vel(ibody)%node(iv)%pos(:)
				call vel_exact(ibody,sphr,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,tmp_vec1)			
				write (*,"(3D15.7)") sphr(1),sphr(2)*180.0/pi,sphr(3)*180.0/pi
				write (*,"(3D15.7)") rhs_check(np+iv), rhs_check(np+nv+iv), rhs_check(np+2*nv+iv)
				write (*,"(3D15.7)") tmp_vec1
				read (*,*) 
			enddo
#endif
			deallocate(coeff,coeff_tmp)
		enddo
		deallocate(p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm)
	end subroutine exact_field


	subroutine vel_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,vec)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: vec(:)

		call ur_exact    (ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,vec(1))
		call utheta_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,vec(2))
		call uphi_exact  (ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,vec(3))
	end subroutine vel_exact


	! P AND DERIVATIVES ^^^^^^^^^^^^^^^^^^^^^^
	subroutine p_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: pres

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		pres = 0d0
		do n=1, nc
			tmp1 = (s**n - n*(2d0*n-1d0)/(2d0*(n+1d0))*s**(-n-1d0)) * (rhof*vis**2/(radbdy(ibody)*dy)**2)
			tmp2 = (-n*(4d0*n**2-1d0)/(n+1d0)*s**(-n-1d0))          * (rhof*vis**2/(radbdy(ibody)*dy)**2)
			! Pnm & P~nm
			do m=0, n
				pres = pres + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo

			! PHInm & PHI~nm
			do m=0, n
				pres = pres + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo
		enddo
	end subroutine p_exact

	subroutine dpdr_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: pres

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		pres = 0d0
		do n=1, nc
			tmp1 = (n* s**(n-1d0) - n*(2d0*n-1d0)/(2d0*(n+1d0))* (-n-1d0)* s**(-n-2d0)) * (rhof*vis**2/(radbdy(ibody)*dy)**2) /(radbdy(ibody)*dy)
			tmp2 = (-n*(4d0*n**2-1d0)/(n+1d0)* (-n-1d0)* s**(-n-2d0))                   * (rhof*vis**2/(radbdy(ibody)*dy)**2) /(radbdy(ibody)*dy)
			! Pnm & P~nm
			do m=0, n
				pres = pres + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo

			! PHInm & PHI~nm
			do m=0, n
				pres = pres + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo
		enddo
	end subroutine dpdr_exact

	subroutine dp_rdtheta_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: pres

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		pres = 0d0
		do n=1, nc
			tmp1 = (s**(n-1d0) - n*(2d0*n-1d0)/(2d0*(n+1d0))*s**(-n-2d0)) * (rhof*vis**2/(radbdy(ibody)*dy)**2) / (radbdy(ibody)*dy)
			tmp2 = (-n*(4d0*n**2-1d0)/(n+1d0)*s**(-n-2d0))                * (rhof*vis**2/(radbdy(ibody)*dy)**2) / (radbdy(ibody)*dy)
			! Pnm & P~nm
			do m=0, n
				pres = pres + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				pres = pres + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo
		enddo
	end subroutine dp_rdtheta_exact

	subroutine dp_rsindphi_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,pres)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: pres

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		pres = 0d0
		do n=1, nc
			tmp1 = (s**(n-1d0) - n*(2d0*n-1d0)/(2d0*(n+1d0))*s**(-n-2d0)) * (rhof*vis**2/(radbdy(ibody)*dy)**2) / (radbdy(ibody)*dy)
			tmp2 = (-n*(4d0*n**2-1d0)/(n+1d0)*s**(-n-2d0))                * (rhof*vis**2/(radbdy(ibody)*dy)**2) / (radbdy(ibody)*dy)
			! Pnm & P~nm
			do m=0, n
				pres = pres + tmp1 * ((-m)*p_nm(n,m)*dsin(m*phi) + m*pt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta))/dsin(theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				pres = pres + tmp2 * ((-m)*f_nm(n,m)*dsin(m*phi) + m*ft_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta))/dsin(theta)
			enddo
		enddo
	end subroutine dp_rsindphi_exact
	! -----------------------------------------------

	! U_r AND DERIVATIVES ^^^^^^^^^^^^^^^^^^^^^^
	subroutine ur_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,ur)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: ur

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		ur = 0d0
		do n=1, nc
			tmp1 = -(n/2d0/(2d0*n+3d0)*s**(n+1d0)- n/4d0*s**(-n) + n*(2d0*n+1d0)/4d0/(2d0*n+3d0)*s**(-n-2d0)) * vis/(rhof*radbdy(ibody)*dy)
			tmp2 = -(n*s**(n-1d0) - n*(2d0*n+1d0)/2d0*s**(-n) + n*(2d0*n-1d0)/2d0*s**(-n-2d0))                * vis/(rhof*radbdy(ibody)*dy)
			! Pnm & P~nm
			do m=0, n
				ur = ur + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo

			! PHInm & PHI~nm
			do m=0, n
				ur = ur + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo
		enddo
	end subroutine ur_exact

	subroutine d2urdr2_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,ur)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: ur

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		ur = 0d0
		do n=1, nc
			tmp1 =-(n/2d0/(2d0*n+3d0)* (n+1d0)*(n+2d0)*s**(n-1d0) - n/4d0*             (-n)*(-n+1d0)*s**(-n-2d0) + n*(2d0*n+1d0)/4d0/(2d0*n+3d0)* (-n-2d0)*(-n-1d0)*s**(-n-4d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			tmp2 =-(n*                 (n-1d0)*(n)    *s**(n-3d0) - n*(2d0*n+1d0)/2d0* (-n)*(-n+1d0)*s**(-n-2d0) + n*(2d0*n-1d0)/2d0*             (-n-2d0)*(-n-1d0)*s**(-n-4d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			! Pnm & P~nm
			do m=0, n
				ur = ur + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo

			! PHInm & PHI~nm
			do m=0, n
				ur = ur + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo
		enddo
	end subroutine d2urdr2_exact
	!------------------------------------------------------------

	! U_theta AND DERIVATIVES ^^^^^^^^^^^^^^^^^^^^^^
	subroutine utheta_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,utheta)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: utheta

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		

		utheta = 0d0
		do n=1, nc
			tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)*s**(n+1d0)+(n-2d0)/4d0/(n+1d0)*s**(-n)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)*s**(-n-2d0)) * vis/(rhof*radbdy(ibody)*dy)
			tmp2 = (s**(n-1d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*s**(-n)-n*(2d0*n-1d0)/2d0/(n+1d0)*s**(-n-2d0)) * vis/(rhof*radbdy(ibody)*dy)
			tmp3 = (s**n-s**(-n-1d0)) * vis/(rhof*radbdy(ibody)*dy)
			! Pnm & P~nm
			do m=0, n
				utheta = utheta + tmp1 * (p_nm(n,m)*dcos(m*phi)+pt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				utheta = utheta + tmp2 * (f_nm(n,m)*dcos(m*phi)+ft_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				utheta = utheta + tmp3 * ((-m)*x_nm(n,m)*dsin(m*phi)+m*xt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo
		enddo
	end subroutine utheta_exact



	subroutine rdutheta_rdr_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,utheta)
		! r d(u_theta/r)/dr
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: utheta

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		

		utheta = 0d0
		do n=1, nc
			tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)*     (n)*s**(n)    +(n-2d0)/4d0/(n+1d0)*                 (-n-1d0)*s**(-n-1d0)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)* (-n-2d0)*s**(-n-3d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)
			tmp2 = (                                 (n-2d0)*s**(n-2d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*     (-n-1d0)*s**(-n-1d0)-n*(2d0*n-1d0)/2d0/(n+1d0)*             (-n-2d0)*s**(-n-3d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)
			tmp3 = (                                 (n-1d0)*s**(n-1d0)-                                     (-n-2d0)*s**(-n-2d0))                                                             * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)
			! Pnm & P~nm
			do m=0, n
				utheta = utheta + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				utheta = utheta + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				utheta = utheta + tmp3 * ((-m)*x_nm(n,m)*dsin(m*phi) + m*xt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo
		enddo
	end subroutine rdutheta_rdr_exact



	subroutine d2uthetadr2_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,utheta)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: utheta

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		utheta = 0d0
		do n=1, nc
			tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)* (n+1d0)*(n+2d0)*s**(n-1d0)    +(n-2d0)/4d0/(n+1d0)*            (-n)*(-n+1d0)*s**(-n-2d0)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)* (-n-2d0)*(-n-1d0)*s**(-n-4d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			tmp2 = (                                     (n-1d0)*(n)*s**(n-3d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*    (-n)*(-n+1d0)*s**(-n-2d0)-n*(2d0*n-1d0)/2d0/(n+1d0)*             (-n-2d0)*(-n-1d0)*s**(-n-4d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			tmp3 = (                                       n*(n+1d0)*s**(n-2d0)-                                    (-n-1d0)*(-n)*s**(-n-3d0))                                             * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			! Pnm & P~nm
			do m=0, n
				utheta = utheta + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				utheta = utheta + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				utheta = utheta + tmp3 * ((-m)*x_nm(n,m)*dsin(m*phi) + m*xt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo
		enddo
	end subroutine d2uthetadr2_exact
	!--------------------------------------------------------

	! U_phi AND DERIVATIVES ^^^^^^^^^^^^^^^^^^^^^^
	subroutine uphi_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,uphi)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: uphi

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		

		uphi = 0d0
		do n=1, nc
			tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)*s**(n+1d0)+(n-2d0)/4d0/(n+1d0)*s**(-n)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)*s**(-n-2d0)) &
				& * vis/(rhof*radbdy(ibody)*dy)
			tmp2 = (s**(n-1d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*s**(-n)-n*(2d0*n-1d0)/2d0/(n+1d0)*s**(-n-2d0)) * vis/(rhof*radbdy(ibody)*dy)
			tmp3 =-(s**n-s**(-n-1d0)) * vis/(rhof*radbdy(ibody)*dy)

			! Pnm & P~nm
			do m=0, n
				uphi = uphi + tmp1 * ((-m)*p_nm(n,m)*dsin(m*phi) + m*pt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				uphi = uphi + tmp2 * ((-m)*f_nm(n,m)*dsin(m*phi) + m*ft_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				uphi = uphi + tmp3 * (x_nm(n,m)*dcos(m*phi) + xt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo
		enddo
	end subroutine uphi_exact

	subroutine d2uphidr2_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,uphi)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: uphi

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		

		uphi = 0d0
		do n=1, nc
			tmp1 = ((n+3d0)/2d0/(n+1d0)/(2d0*n+3d0)* (n+1d0)*(n+2d0)*s**(n-1d0)    +(n-2d0)/4d0/(n+1d0)*            (-n)*(-n+1d0)*s**(-n-2d0)-n*(2d0*n+1d0)/4d0/(n+1d0)/(2d0*n+3d0)* (-n-2d0)*(-n-1d0)*s**(-n-4d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			tmp2 = (                                     (n-1d0)*(n)*s**(n-3d0)+(n-2d0)*(2d0*n+1d0)/2d0/(n+1d0)*    (-n)*(-n+1d0)*s**(-n-2d0)-n*(2d0*n-1d0)/2d0/(n+1d0)*             (-n-2d0)*(-n-1d0)*s**(-n-4d0)) * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			tmp3 =-(                                       n*(n+1d0)*s**(n-2d0)-                                    (-n-1d0)*(-n)*s**(-n-3d0))                                             * vis/(rhof*radbdy(ibody)*dy) / (radbdy(ibody)*dy)**2d0
			! Pnm & P~nm
			do m=0, n
				uphi = uphi + tmp1 * ((-m)*p_nm(n,m)*dsin(m*phi) + m*pt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo

			! PHInm & PHI~nm
			do m=0, n
				uphi = uphi + tmp2 * ((-m)*f_nm(n,m)*dsin(m*phi) + m*ft_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta)) / dsin(theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				uphi = uphi + tmp3 * (x_nm(n,m)*dcos(m*phi) + xt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo
		enddo
	end subroutine d2uphidr2_exact
	!---------------------------------------------

	subroutine omega_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega_vec)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: omega_vec(:)

		call omegar_exact    (ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega_vec(1))
		call omegatheta_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega_vec(2))
		call omegaphi_exact  (ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega_vec(3))
	end subroutine omega_exact

	! Omega_r ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	subroutine omegar_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: omega

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		omega = 0d0
		do n=1, nc
			tmp1 = n*(n+1d0)*(s**(n-1d0)-s**(-n-2d0)) * vis/(rhof*(radbdy(ibody)*dy)**2)
			! Xnm & X~nm
			do m=0, n
				omega = omega + tmp1 * (x_nm(n,m)*dcos(m*phi) + xt_nm(n,m)*dsin(m*phi)) * plgndr(n,m,dcos(theta))
			enddo
		enddo
	end subroutine omegar_exact
	!----------------------------------------------

	! Omega_theta ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	subroutine omegatheta_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: omega

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2, tmp3

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		omega = 0d0
		do n=1, nc
			tmp1 = -(s**n/(n+1d0)+(2d0*n-1d0)/2d0/(n+1d0)*s**(-n-1d0)) * vis/(rhof*(radbdy(ibody)*dy)**2)
			tmp2 = -(4d0*n**2d0-1d0)/(n+1d0)*s**(-n-1d0)               * vis/(rhof*(radbdy(ibody)*dy)**2)
			tmp3 = (n+1d0)*s**(n-1d0) + n*s**(-n-2d0)                  * vis/(rhof*(radbdy(ibody)*dy)**2)
			! Pnm & P~nm
			do m=0, n
				omega = omega + tmp1 * (-m*p_nm(n,m)*dsin(m*phi) + m*pt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta))/dsin(theta)
			enddo

			! fnm & f~nm
			do m=0, n
				omega = omega + tmp2 * (-m*f_nm(n,m)*dsin(m*phi) + m*ft_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta))/dsin(theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				omega = omega + tmp3 * (x_nm(n,m)*dcos(m*phi) + xt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo
		enddo
	end subroutine omegatheta_exact
	!----------------------------------------------

	! Omega_phi ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	subroutine omegaphi_exact(ibody,coord,p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm,nc,omega)
		implicit none
		integer, intent(in) :: ibody, nc
		real(8), intent(in) :: coord(ndim)
		real(8), intent(in), dimension(1:nc,0:nc) :: p_nm,pt_nm,f_nm,ft_nm,x_nm,xt_nm
		real(8), intent(out) :: omega

		integer :: n, m
		real(8) :: rad, theta, phi
		real(8) :: tmp1, tmp2, tmp3

		rad = coord(1)
		theta = coord(2)
		phi = coord(3)

		s = coord(1)/radbdy(ibody)
		
		omega = 0d0
		do n=1, nc
			tmp1 = (s**n/(n+1d0)+(2d0*n-1d0)/2d0/(n+1d0)*s**(-n-1d0)) * vis/(rhof*(radbdy(ibody)*dy)**2)
			tmp2 = (4d0*n**2d0-1d0)/(n+1d0)*s**(-n-1d0)               * vis/(rhof*(radbdy(ibody)*dy)**2)
			tmp3 = (n+1d0)*s**(n-1d0) + n*s**(-n-2d0)                 * vis/(rhof*(radbdy(ibody)*dy)**2)
			! Pnm & P~nm
			do m=0, n
				omega = omega + tmp1 * (p_nm(n,m)*dcos(m*phi) + pt_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! fnm & f~nm
			do m=0, n
				omega = omega + tmp2 * (f_nm(n,m)*dcos(m*phi) + ft_nm(n,m)*dsin(m*phi)) * dplgndr(n,m,theta)
			enddo

			! Xnm & X~nm
			do m=0, n
				omega = omega + tmp3 * (-m*x_nm(n,m)*dsin(m*phi) + m*xt_nm(n,m)*dcos(m*phi)) * plgndr(n,m,dcos(theta))/dsin(theta)
			enddo
		enddo
	end subroutine omegaphi_exact
	!----------------------------------------------

	subroutine form_cage_nodes
		implicit none

		write (*,*) "FORMING THE CAGE NODES"
		allocate(cage_pres(nbody), cage_vel(nbody))
		cage_pres%num = 0
		cage_vel%num  = 0

		call initialize_gridvertindex (nx, my, mz)

		do m=1, nbody
			allocate(tmp_ind(int(2*(2*radbdy(m))**3)))

			! FINDING PRESSURE CAGE NODES
			tmp_ind = 0
			tmp_num = 0
			do ijk = 1, nvert
				call ind1t3(ijk,i,j,k)
				call cell_pres_dist(i,j,k,m,tmp_vec1,dist1)

				if (radbdy(m)<dist1.and.dist1<=radbdy(m)+2*dr) then	! FINDING THE NODE BEING NEAR THE PARTICLE SURFACE IN FLUID PHASE
					im = i-1
					ip = i+1
					jm = j-1
					jp = j+1
					km = k-1
					kp = k+1

					! CHECKING IF THE FLUID PRESSURE NODE HAS A NEIGHBOR IN SOLID
					neighb_insolid = .false.
					do kk=km, kp
						do jj=jm, jp
							do ii=im, ip
								iii = ii
								jjj = jj
								kkk = kk
#if !PARALLEL
								if (iii<1)  iii=iii+nx
								if (iii>nx) iii=iii-nx
#endif
								if (jjj<1)  jjj=jjj+my
								if (jjj>my) jjj=jjj-my

								if (kkk<1)  kkk=kkk+mz
								if (kkk>mz) kkk=kkk-mz

								call cell_pres_dist(iii,jjj,kkk,m,tmp_vec1,dist2)

								if (dist2<=radbdy(m)) neighb_insolid = .true.

								if (neighb_insolid) exit
							enddo
							if (neighb_insolid) exit
						enddo
						if (neighb_insolid) exit
					enddo
					! ADDING THE NODE TO THE LIST
					if (neighb_insolid) then
						! NEGLECTING THE POINTS WHICH HAVE THETA = 0.0
						call cell_pres_dist(i,j,k,m,tmp_vec1,dist1)
						cart(:) = tmp_vec1(:)
						call cart_to_sphr(cart,m,sphr)

						if (.not.((abs(sphr(2))<1d-3).or.(abs(sphr(2)-pi*half)<1d-3))) then
							tmp_num = tmp_num+1
							tmp_ind(tmp_num) = ijk
						endif
!						tmp_num = tmp_num+1
!						tmp_ind(tmp_num) = ijk
					endif
				endif
			enddo

			cage_pres(m)%num = tmp_num
			allocate(cage_pres(m)%node(tmp_num))

			do n=1, tmp_num
				call ind1t3(tmp_ind(n),i,j,k)
				call cell_pres_dist(i,j,k,m,tmp_vec1,dist1)
				cage_pres(m)%node(n)%ind = tmp_ind(n)
				cage_pres(m)%node(n)%pos(:) = tmp_vec1(:)
			enddo

			! EXCLUDING THE NODES WHICH ARE CLOSER TO PARTICLE CENTER COMPARE TO THE PRESSURE POSITION
			tmp_num = 0
			tmp_ind = 0
			do n=1, cage_pres(m)%num
				call ind1t3(cage_pres(m)%node(n)%ind,i,j,k)
				call cell_pres_dist(i,j,k,m,tmp_vec1,dist1)

				ip = i+1
				jp = j+1
				kp = k+1

				! CHECKING FOR ALL VERTICES OFF THE CELL
				do ii=i, ip
					do jj=j, jp
						do kk=k, kp
							iii = ii
							jjj = jj
							kkk = kk
#if !PARALLEL
							if (iii<1)  iii=iii+nx
							if (iii>nx) iii=iii-nx
#endif
							if (jjj<1)  jjj=jjj+my
							if (jjj>my) jjj=jjj-my

							if (kkk<1)  kkk=kkk+mz
							if (kkk>mz) kkk=kkk-mz

							call cell_vertex_dist(iii,jjj,kkk,m,tmp_vec1,dist2)
							if (dist2<=dist1) then
								call vertindex_ijk(ijk,iii,jjj,kkk,nx,my,mz)

								registered = .false.
								do o=1, tmp_num
									if (tmp_ind(o)==ijk) then
										registered = .true.
										exit
									endif
								enddo
								if (.not.registered) then
									tmp_num = tmp_num + 1
									tmp_ind(tmp_num) = ijk
								endif
							endif
						enddo
					enddo
				enddo
			enddo
			exclude_num = tmp_num
			allocate(exclude_ind(tmp_num))
			exclude_ind(1:tmp_num) = tmp_ind(1:tmp_num)


			! FINDING VELOCITY CAGE NODES
			tmp_ind = 0
			tmp_num = 0
			do n=1, cage_pres(m)%num
				! TO FIGURE OUT THE CORRESPONDIG VELOCITY CELL, THE PRESSURE CELL POSITION IS EXTENDED IN R-COORDINATE
				! A LITTLE AND THE POSITION OF NEW POINT DEFINES THE MOST SUITABLE CELL TO LOOK FOR VELOCITY NODES
				call ind1t3(cage_pres(m)%node(n)%ind,i,j,k)
				cart = cage_pres(m)%node(n)%pos
				call cart_to_sphr(cart,m,sphr)
				sphr(1) = sphr(1) + 1d-3
				call sphr_to_cart(sphr,m,cart)
				call find_cell(cart,ijk)
				call ind1t3(ijk,i,j,k)

				ip = i+1
				jp = j+1
				kp = k+1

				! LOOKING FOR FARTHUR VERTICES WITH RESPECT TO PRESSURE POSITION
				do ii=i, ip
					do jj=j, jp
						do kk=k, kp
							iii = ii
							jjj = jj
							kkk = kk
#if !PARALLEL
							if (iii<1)  iii=iii+nx
							if (iii>nx) iii=iii-nx
#endif
							if (jjj<1)  jjj=jjj+my
							if (jjj>my) jjj=jjj-my

							if (kkk<1)  kkk=kkk+mz
							if (kkk>mz) kkk=kkk-mz

							call vertindex_ijk(ijk_1,iii,jjj,kkk,nx,my,mz)
							! LOOKING IN EXCLUSIONS LIST
							registered = .false.
							do o=1, exclude_num
								if (exclude_ind(o)==ijk_1) then
									registered = .true.
									exit
								endif
							enddo
	
							if (.not.registered) then
								registered = .false.
								do o=1, tmp_num
									if (tmp_ind(o)==ijk_1) then
										registered = .true.
										exit
									endif
								enddo
								if (.not.registered) then
									! NEGLECTING THE POINTS WHICH HAVE THETA = 0.0
									call cell_pres_dist(iii,jjj,kkk,m,tmp_vec1,dist1)
									cart(:) = tmp_vec1(:)
									call cart_to_sphr(cart,m,sphr)

									if (.not.((abs(sphr(2))<1d-3).or.(abs(sphr(2)-pi*half)<1d-3))) then
										tmp_num = tmp_num + 1
										tmp_ind(tmp_num) = ijk_1
									endif
!									tmp_num = tmp_num + 1
!									tmp_ind(tmp_num) = ijk_1
								endif
							endif
						enddo
					enddo
				enddo
			enddo
			cage_vel(m)%num = tmp_num
			allocate(cage_vel(m)%node(tmp_num))
			do n=1, tmp_num
				call ind1t3(tmp_ind(n),i,j,k)
				call cell_vertex_dist(i,j,k,m,tmp_vec1,dist1)
				cage_vel(m)%node(n)%ind = tmp_ind(n)
				cage_vel(m)%node(n)%pos(:) = tmp_vec1(:)
			enddo
			deallocate(tmp_ind,exclude_ind)
		enddo
	end subroutine form_cage_nodes

	subroutine change_coordinates
		implicit none

		write (*,*) "CHANGING THE COORDINATES TO SPHERICAL"
		do m=1, nbody
			do n=1, cage_pres(m)%num
				cart(:) = cage_pres(m)%node(n)%pos(:)
				call cart_to_sphr(cart,m,sphr)
				cage_pres(m)%node(n)%pos(:) = sphr(:)
			enddo

			do n=1, cage_vel(m)%num
				cart(:) = cage_vel(m)%node(n)%pos(:)
				call cart_to_sphr(cart,m,sphr)
				cage_vel(m)%node(n)%pos(:) = sphr(:)
			enddo
		enddo
	end subroutine change_coordinates

	subroutine cage_output
		implicit none

		write (*,*) "TOTAL NODES = ", nx*my*mz
		write (*,*) "PRESS NODES = ", cage_pres(1)%num
		write (*,*) "VEL   NODES = ", cage_vel(1)%num
		write (*,*) "EXCL. NODES = ", exclude_num

		write (*,*) "WRITING CAGE NODES OUTPUT"

		open (unit=1, file=trim(run_name)//"_cage_pres_N4_post.dat", status="replace")
	!	write (1,*) "variables=x,y"
		do j=1, my
			registered = .false.
			write (1,*)'ZONE T= "', j, ' " '
			do i=1, cage_pres(1)%num
				if (int(cage_pres(1)%node(i)%pos(2)+xc(1,2))==j) then
					write (1,"(2D15.7)") cage_pres(1)%node(i)%pos(1), cage_pres(1)%node(i)%pos(3)
					registered = .true.
				endif
			enddo
			if (.not.registered) write (1,"(2D15.7)") 0d0, 0d0
		enddo
		close (1)

		open (unit=1, file=trim(run_name)//"_cage_vel_N4_post.dat", status="replace")
	!	write (1,*) "variables=x,y"
		do j=1, my
			registered = .false.
			write (1,*)'ZONE T= "', j, ' " '
			do i=1, cage_vel(1)%num
				if (int(cage_vel(1)%node(i)%pos(2)+xc(1,2))==j) then
					write (1,"(2D15.7)") cage_vel(1)%node(i)%pos(1), cage_vel(1)%node(i)%pos(3)
					registered = .true.
				endif
			enddo
			if (.not.registered) write (1,"(2D15.7)") 0d0, 0d0
		enddo
		close (1)

		open (unit=1, file=trim(run_name)//"_grid_N4_post.dat", status="replace")
	!	write (1,*) "variables=x,y"
		do j=1, my
			write (1,*)'ZONE T= "', j, ' " ', "I = ", nx, " J = ", my
			do k=1, mz
				do i=1, nx
					write (1,"(2D15.7)") i-xc(1,1),k-xc(1,3)
				enddo
			enddo
		enddo
		close (1)

		open (unit=1, file=trim(run_name)//"_cage_3D_N4_post.dat", status="replace")
		write (1,*) "variables=x,y,z,surf"
		write (1,*)'ZONE ', "I = ", nx, " J = ", my, " K = ", mz, "F= point"
		do k=1, mz
			do j=1, my
				do i=1, nx
					write (1,"(4D15.7)") i-xc(1,1),j-xc(1,2),k-xc(1,3),10d0
				enddo
			enddo
		enddo
		write (1,*)'ZONE '
		do i=1, cage_pres(1)%num
			write (1,"(4D15.7)") cage_pres(1)%node(i)%pos(1), cage_pres(1)%node(i)%pos(2), cage_pres(1)%node(i)%pos(3), 0d0
		enddo
		write (1,*)'ZONE '
		do i=1, cage_vel(1)%num
			write (1,"(4D15.7)") cage_vel(1)%node(i)%pos(1), cage_vel(1)%node(i)%pos(2), cage_vel(1)%node(i)%pos(3), 1d0
		enddo
		close (1)
	end subroutine cage_output

	subroutine cell_center_dist(i,j,k,m,center,dist)
		implicit none
		integer, intent(in) :: i,j,k,m
		real(8), intent(out) :: center(ndim), dist
		
		center(1) = dble(i) + half
		center(2) = dble(j) + half
		center(3) = dble(k) + half

		center(:) = center(:) - xc(m,:)

		dist = sqrt(dot_product(center,center))
	end subroutine cell_center_dist

	subroutine cell_pres_dist(i,j,k,m,center,dist)
		implicit none
		integer, intent(in) :: i,j,k,m
		real(8), intent(out) :: center(ndim), dist
		
		center(1) = dble(i) + half
		center(2) = dble(j)
		center(3) = dble(k)

		center(:) = center(:) - xc(m,:)
		dist = sqrt(dot_product(center,center))
	end subroutine cell_pres_dist

	subroutine cell_vertex_dist(i,j,k,m,center,dist)
		implicit none
		integer, intent(in) :: i,j,k,m
		real(8), intent(out) :: center(ndim), dist
		
		center(1) = dble(i)
		center(2) = dble(j)
		center(3) = dble(k)

		center(:) = center(:) - xc(m,:)
		dist = sqrt(dot_product(center,center))
	end subroutine cell_vertex_dist

	subroutine cart_to_sphr(cart,m,sphr)
		implicit none
		integer, intent(in) :: m
		real(8), intent(in) :: cart(ndim)
		real(8), intent(out) :: sphr(ndim)
		real(8) :: r, phi, theta, r_tmp

		integer :: idim

!		r = sqrt(dot_product(cart,cart))
!		if (abs(cart(3))<r) then
!			theta = acos(cart(3)/r)
!		else
!			theta = sign(cart(3),cart(3))*acos(1d0)
!		endif
!		r_tmp = r * sin(theta)
!		if (r_tmp>small_number) then
!			if (abs(cart(1))<r_tmp) then
!				phi = acos(abs(cart(1))/r_tmp)
!			else
!				phi = acos(1d0)
!			endif
!
!			if (cart(1)<0d0.and.cart(2)>=0d0) then
!				phi = pi-phi
!			elseif (cart(1)<0d0.and.cart(2)<0d0) then
!				phi = pi+phi
!			elseif (cart(1)>=0d0.and.cart(2)<0d0) then
!				phi = twopi-phi
!			endif
!		else
!			phi = zero
!		endif
!		sphr(1) = r
!		sphr(2) = theta
!		sphr(3) = phi

		r = sqrt(dot_product(cart,cart))
		if (abs(cart(1))<r) then
			theta = acos(-cart(1)/r)
		else
			theta = sign(cart(1),-cart(1))*acos(1d0)
		endif
		r_tmp = r * sin(theta)
		if (r_tmp>small_number) then
			if (abs(cart(2))<r_tmp) then
				phi = acos(abs(cart(2))/r_tmp)
			else
				phi = acos(1d0)
			endif

			if (cart(2)>=0d0.and.cart(3)>=0d0) then

			elseif (cart(2)<0d0.and.cart(3)>=0d0) then
				phi = pi-phi
			elseif (cart(2)<0d0.and.cart(3)<0d0) then
				phi = pi+phi
			elseif (cart(2)>=0d0.and.cart(3)<0d0) then
				phi = twopi-phi
			endif
		else
			phi = zero
		endif
		sphr(1) = r
		sphr(2) = theta
		sphr(3) = phi
	end subroutine cart_to_sphr

	subroutine sphr_to_cart(sphr,m,cart)
		implicit none
		integer, intent(in) :: m
		real(8), intent(in) :: sphr(ndim)
		real(8), intent(out) :: cart(ndim)
		real(8) :: r, phi, theta, r_tmp

		integer :: idim
		
		r     = sphr(1)
		theta = sphr(2)
		phi   = sphr(3)

		r_tmp = r * sin(theta)

!		cart(3) = r     * cos(theta) + xc(m,3)
!		cart(1) = r_tmp * cos(phi)   + xc(m,1)
!		cart(2) = r_tmp * sin(phi)   + xc(m,2)

		cart(1) = -r     * cos(theta) + xc(m,1)
		cart(2) =  r_tmp * cos(phi)   + xc(m,2)
		cart(3) =  r_tmp * sin(phi)   + xc(m,3)
	end subroutine sphr_to_cart

	subroutine vec_cart_to_sphr(vec,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: vec(ndim)

		real(8) :: tmp(3)
		real(8) :: rotate(3,3)

		integer :: dim1,dim2,dim3,dim4

		!ROTATION OF THE STRESS TENSOR ALIGNED WITH THE SPHERICAL COORDINATE
		rotate(1,1) =-cos(th)
		rotate(2,1) = sin(th) * cos(phi)
		rotate(3,1) = sin(th) * sin(phi)

		rotate(1,2) = sin(th)
		rotate(2,2) = cos(th) * cos(phi)
		rotate(3,2) = cos(th) * sin(phi)

		rotate(1,3) = 0d0
		rotate(2,3) =-sin(phi)
		rotate(3,3) = cos(phi)

		tmp = 0d0
		do dim1=1, ndim
			do dim2=1, ndim
				tmp(dim1) = tmp(dim1) + rotate(dim2,dim1) * vec(dim2)
			enddo
		enddo
		vec(:) = tmp(:)
	end subroutine vec_cart_to_sphr

	subroutine vec_sphr_to_cat(vec,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: vec(ndim)

		real(8) :: tmp(3)
		real(8) :: rotate(3,3)

		integer :: dim1,dim2,dim3,dim4

		!ROTATION OF THE STRESS TENSOR ALIGNED WITH THE SPHERICAL COORDINATE
		rotate(1,1) =-cos(th)
		rotate(2,1) = sin(th) * cos(phi)
		rotate(3,1) = sin(th) * sin(phi)

		rotate(1,2) = sin(th)
		rotate(2,2) = cos(th) * cos(phi)
		rotate(3,2) = cos(th) * sin(phi)

		rotate(1,3) = 0d0
		rotate(2,3) =-sin(phi)
		rotate(3,3) = cos(phi)

		tmp = 0d0
		do dim1=1, ndim
			do dim2=1, ndim
				tmp(dim1) = tmp(dim1) + rotate(dim1,dim2) * vec(dim2)
			enddo
		enddo
		vec(:) = tmp(:)
	end subroutine vec_sphr_to_cat

	subroutine find_cell(vec,ind)
		implicit none
		real(8), intent(in)  :: vec(ndim)
		integer, intent(out) :: ind
		integer :: i, j, k, ijk

		i = int(vec(1))
		j = int(vec(2))
		k = int(vec(3))

		call vertindex_ijk(ijk,i,j,k,nx,my,mz)
		ind = ijk
	end subroutine find_cell
end subroutine physalis

end module physalis_mod


SUBROUTINE  xwnnls(vprime,aprime,gamma,mv, nv)
	use global_data, only : prcn
	IMPLICIT NONE
!!$    This subroutine is used to estimate gamma by casting A' = -
    !!\gamma v' into a  linear least squares problem. The approximate
    !! set of equations to be solved are                             
    !!      V'\Gamma = -A' 
!!$  V' => mv X nv,  \Gamma => nv X 1,   A' => mv X 1
!!$ The equations are solved with the non negativity constraints in
    !! the last few rows of the \Gamma vector.

	INTEGER , INTENT(in) :: mv,nv
	INTEGER :: me,lw,liw,L,mdw,mode, row
	REAL(prcn), Intent(out) :: gamma(nv)
	REAL(prcn) :: prgopt(1),upval
	REAL(prcn),INTENT(inout) ::  vprime(mv,nv), aprime(mv)
	REAL(prcn)  :: rnorm
	REAL(prcn) :: vprime_temp(mv,nv), aprime_temp(mv),gamma_temp(nv)
	INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
	REAL(prcn), DIMENSION(:,:),ALLOCATABLE ::  E
	REAL(prcn), DIMENSION(:),ALLOCATABLE ::  F
	REAL(prcn), DIMENSION(:), ALLOCATABLE :: work
	REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: W

	me = 0 ! For the input to DWNNLS

!!$    Interchange rows and columns corresponding to \gamma_11, 
    !!\gamma_22 and \gamma_33 (since non negativity constraints are
    !! imposed on these values)

	vprime_temp = vprime
!!$    do row = 1, mv
!!$       vprime(row,1) = vprime_temp(row,nv-ndim+1) ! v'_i1 --> v'_i7
!!$       vprime(row,nv-ndim+1) = vprime_temp(row,1) ! v'_i7 --> v'_i1
!!$       vprime(row,1+ndim+1) = vprime_temp(row,nv-ndim+2) ! v'_i5 --> v'_i8
!!$       vprime(row, nv-ndim+2) = vprime_temp(row,1+ndim+1) ! v'_i8 --> v'_i5
!!$    end do
	aprime_temp = aprime
	mdw  = me+mv
	!k = max(ma+mg,nmat)
	lw = me+mv+5*nv
	liw = me+mv+nv

	IF(.NOT.ALLOCATED(iwork)) ALLOCATE(iwork(liw))
	IF(.NOT.ALLOCATED(work)) ALLOCATE(work(mv+5*nv))
	IF(.NOT.ALLOCATED(W)) ALLOCATE(W(mdw,nv+1))
	IF(me.GT.0) THEN 
		W(1:me,1:nv) = E(1:me,1:nv)
		W(1:me,nv+1) = F(1:me)
	ENDIF
	W(me+1:me+mv,1:nv) = vprime(1:mv,1:nv)
	W(me+1:me+mv,nv+1) = aprime(1:mv)
	iwork(1) = lw
	iwork(2) = liw
	prgopt(1) = 1.
	!PRINT*,'MA = ', MV
	!PRINT*,'NA = ', NV
	!L = 3
	!L=nv-ndim
	L = nv
	!PRINT*,'L = ', L
	CALL DWNNLS (W, MDW, ME, MV, NV, L, PRGOPT, gamma_temp, RNORM, MODE,IWORK, WORK)
	gamma(1:nv) = gamma_temp(1:nv)
	!IF(mode.NE.0)   
	!PRINT*,'Mode=',mode

!!$    gamma(1) = gamma_temp(nv-ndim+1)
!!$    gamma(nv-ndim+1) = gamma_temp(1)
!!$    gamma(1+ndim+1) = gamma_temp(nv-ndim+2)
!!$    gamma(nv-ndim+2) = gamma_temp(1+ndim+1)
!	do row=1,nv
!		PRINT*,'gamma = ' ,gamma(row)
!	end do

  end SUBROUTINE xwnnls





	
#if 0
	subroutine calc_ugrad
		implicit none

		do i=1, nx
			if (I_AM_NODE_ZERO.and.mod(i,10)==0) write (*,*) " PLANE # = ", i
			ip=i+1
			im=i-1
#if !PARALLEL
			if (ip>nx) ip=1
			if (im<1)  im=nx
#endif
			do dir=1, ndim
				do idim=1, ndim
					do k=1, mz
						do j=1, my2
							if (dir==1) then
								uf1(j,k) = (u(ip,j,k,idim)-u(im,j,k,idim)) / (two*dx)
							elseif (dir==2) then
								uf1(j,k) = u(i,j,k,idim) * wy(j)
							elseif (dir==3) then
								uf1(j,k) = u(i,j,k,idim) * wz(k)
							endif
						enddo
					enddo
					call ff2cr(uf1,ur1)
					ugrad(i,:,:,idim,dir) = ur1(:,:)
				enddo
			enddo
		enddo

	
!		open (unit=1, file=trim(run_name)//"_ugrad_N4_post.dat", status="replace")
!		write (1,*) "variables=i,j,k,u11,u12,u13" !,u21,u22,u23,u31,u32,u33"
!		write (1,*) "zone  i=",nx," j=",my," k=",mz," f=point"
!
!		do k=1,mz		
!			do j=1, my
!				do i=1, nx
!!					do idim=1, ndim
!						write (1,"(3I6,3D15.7)") i,j,k,ugrad(i,j,k,1,:)
!!					enddo
!				enddo
!			enddo
!		enddo
!		close(1)


	end subroutine calc_ugrad

	subroutine calc_eps
		implicit none
		integer :: iii, jjj, kkk

!		do dim3=1, ndim
			do dim2=1, ndim
				do dim1=1, ndim
					do k=1, mz
						do j=1, my
							do i=1, nx
								if (fluid_atijk(i,j,k)) then
!									epsij_loc(dim1,dim2) = epsij_loc(dim1,dim2) + 2 * vis * ugrad(i,j,k,dim1,dim3) * ugrad(i,j,k,dim2,dim3)
									
									tmp1 = vis * (ubcp(i,j,k,dim1)-ufmean(dim1)) * diffn(i,j,k,dim2)
									tmp2 = vis * (ubcp(i,j,k,dim2)-ufmean(dim2)) * diffn(i,j,k,dim1)
									epsij_loc(dim1,dim2) = epsij_loc(dim1,dim2) + tmp1 + tmp2
	

									im = i-1
									ip = i+1
									jm = i-1
									jp = i+1
									km = i-1
									kp = i+1

									neighb_insolid = .false.
									do kk=km, kp
										do jj=jm, jp
											do ii=im, ip
												iii = ii
												jjj = jj
												kkk = kk
#if !PARALLEL
												if (iii<1)  iii=nx
												if (iii>nx) iii=1
#endif
												if (jjj<1)  jjj=my
												if (jjj>my) jjj=1
												if (kkk<1)  kkk=my
												if (kkk>mz) kkk=1

												if (.not.fluid_atijk(iii,jjj,kkk)) neighb_insolid = .true.
												if (neighb_insolid) exit
											enddo
											if (neighb_insolid) exit
										enddo
										if (neighb_insolid) exit
									enddo

!									if (.not.neighb_insolid) epsij_far_loc(dim1,dim2) = epsij_far_loc(dim1,dim2) + 2 * vis * ugrad(i,j,k,dim1,dim3) * ugrad(i,j,k,dim2,dim3)

									if (.not.neighb_insolid) epsij_far_loc(dim1,dim2) = epsij_far_loc(dim1,dim2) + tmp1 + tmp2
								endif
							enddo
						enddo
					enddo
				enddo
			enddo
!		enddo

#if PARALLEL
	epsij = 0d0
	GLOBAL_DOUBLE_SUM(epsij_loc,epsij,9,decomp_group)

	epsij_far = 0d0
	GLOBAL_DOUBLE_SUM(epsij_far_loc,epsij,9,decomp_group)
#else
	epsij = epsij_loc
	epsij_far = epsij_far_loc
#endif
	epsij = epsij/count_fluid/(tke*umeanslip/dia_phys)
	eps   = (epsij(1,1)+epsij(2,2)+epsij(3,3))/2

	epsij_far = epsij_far/count_fluid/(tke*umeanslip/dia_phys)
	eps_far   = (epsij_far(1,1)+epsij_far(2,2)+epsij_far(3,3))/2

	end subroutine calc_eps


	subroutine transfer_velocity
		implicit none
		integer :: left, right, l
		real(8), allocatable :: trans_s(:), trans_r(:)
#if PARALLEL
		left = myid-1
		right = myid+1

		if (left<0) left=left+nproc
		if (right>nproc-1) right=right-nproc

		allocate(trans_s(my*mz*3), trans_r(my*mz*3))

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					trans_s(l)=ubcp(1,j,k,idim)
				enddo
			enddo
		enddo

	
		call MPI_SENDRECV(trans_s,my*mz*3,MPI_DOUBLE_PRECISION,left,0,trans_r,my*mz*3,MPI_DOUBLE_PRECISION,right,0,decomp_group,status,err_code)

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					ubcp(nx+1,j,k,idim)=trans_r(l)
				enddo
			enddo
		enddo

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					trans_s(l)=ubcp(nx,j,k,idim)
				enddo
			enddo
		enddo
	
		call MPI_SENDRECV(trans_s,my*mz*3,MPI_DOUBLE_PRECISION,right,1,trans_r,my*mz*3,MPI_DOUBLE_PRECISION,left,1,decomp_group,status,err_code)

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					ubcp(0,j,k,idim)=trans_r(l)
				enddo
			enddo
		enddo
#else
		ubcp(nx+1,:,:,:) = ubcp(1,:,:,:)
#endif	
	end subroutine transfer_velocity

#endif

