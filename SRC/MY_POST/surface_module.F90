module surface_module
#include "../FLO/ibm.h"
	use global_data

	use nlmainarrays , only : ubcp, pbcp
	use constants
	use general_funcs
	use string_funcs
!	use boundary_condition
	use dependent_functions
!	use bcsetarrays
	implicit none

contains

	subroutine surface_field(input_name)
		implicit none

		type :: phi_data
			real(8) :: ph
			real(8) :: eps_visc, eps_pres
	!		real(8), dimension(:,:), pointer :: visc_comp
			real(8), dimension(1)    :: pres_comp

			real(8), dimension(ndim) :: pgrad,vel,pos
			real(8), dimension(ndim) :: div2vel
		endtype phi_data

		type :: theta_data
			real(8) :: th
			integer :: nphi
			type (phi_data), dimension(:), allocatable :: phi
		endtype theta_data
		
		type :: surface_stress
			integer :: nth
			type (theta_data), dimension(:), pointer :: theta
		endtype surface_stress

		type(surface_stress) :: in1, surf, out1

		character*50 input_name

		real(prcn), allocatable :: vel_th(:,:,:), pres_th(:,:)


		real(prcn), pointer :: out1_vel(:), surf_vel(:), in1_vel(:)
		real(prcn), pointer :: out1_pres, out2_pres, surf_pres, in1_pres
	!	real(8), dimension(:,:), pointer :: bndarray

		real(prcn) :: tmp1, rad, r, ri, ro, x, y, z, dth, th, dphi, phi, dist
		integer :: ith, nth, thp, thm, iphi, nphi, phim, phip
		integer :: iphs, nbnd, idim
		integer :: i, j, k, l, m, n, dim1, dim2

		integer :: ib, ie, jb, je, kb, ke, onew, ii, jj, kk
	
		integer, dimension(ndim) :: vcellb, vcelli, vcello
		integer, dimension(ndim) :: pcellb, pcelli, pcello
		real(prcn), dimension(ndim) :: xl, xli, xlo, xx
		real(prcn), dimension(ndim) :: xlp, xlip, xlop
		real(prcn), dimension(ndim) :: ul, uli, ulo, ul_p, ulo_p, v_tmp
		real(prcn) :: pl, pii, plo
		real(prcn), dimension(ndim) :: pgl, pgli, pglo

		real(prcn) :: rdummy(ndim), umagnorm, umagnormo, umagnormo2
		real(prcn) :: ddr, ddth, ddphi

		logical :: velterm

		character*50 filename


		nth  = 180
		nphi = 360

		dth = pi/nth
		dphi= twopi/nphi

		allocate(vel_th(3,0:nth,ndim))	
		allocate(pres_th(3,0:nth))
	
		allocate(in1%theta(0:nth))
		allocate(surf%theta(0:nth))
		allocate(out1%theta(0:nth))

		do ith=0, nth
			allocate(in1%theta(ith)%phi(0:nphi-1))
			allocate(surf%theta(ith)%phi(0:nphi-1))
			allocate(out1%theta(ith)%phi(0:nphi-1))
		enddo

		vel_th = 0d0
		pres_th = 0d0

		do m=1, nbody
			ddr  = dr * dx
			r   = radbdy(m) * dx
			ri  = r-ddr
			ro  = r+ddr	! * dx

			do ith=0, nth
				th=ith*dth
				do iphi=0, nphi-1
					phi=iphi*dphi

					xx(1) =-cos(th)
					xx(2) = sin(th)*cos(phi)
					xx(3) = sin(th)*sin(phi)
			
					do n=1, ndim
						xl(n)=xc(m,n)+ xx(n)*radbdy(m)
						surf%theta(ith)%phi(iphi)%pos(n) = xl(n)
						ul(n)=zero

						xli(n)=xc(m,n)+ xx(n)*radibdy(m)
						in1%theta(ith)%phi(iphi)%pos(n) = xli(n)
						uli(n)=zero

						xlo(n)=xc(m,n)+ xx(n)*radobdy(m)
						out1%theta(ith)%phi(iphi)%pos(n) = xx(n)
						ulo(n)=zero
					enddo

					xlp(1)   = xl(1)-0.5
					xlp(2:3) = xl(2:3)

					xlip(1)   = xli(1)-0.5
					xlip(2:3) = xli(2:3)

					xlop(1)   = xlo(1)-0.5
					xlop(2:3) = xlo(2:3)

					do n=1, ndim
						if (xl(n).lt.zero) then 
							vcellb(n) = int(xl(n)-1)
						else 
							vcellb(n) = int(xl(n))
						endif
						if (xlp(n).lt.zero) then 
							pcellb(n) = int(xlp(n)-1)
						else 
							pcellb(n) = int(xlp(n))
						end if

						if (xli(n).lt.zero) then 
							vcelli(n) = int(xli(n)-1)
						else 
							vcelli(n) = int(xli(n))
						endif
						if (xlip(n).lt.zero) then 
							pcelli(n) = int(xlip(n)-1)
						else 
							pcelli(n) = int(xlip(n))
						end if

						if(xlo(n).lt.zero) then 
							vcello(n) = int(xlo(n)-1)
						else 
							vcello(n) = int(xlo(n))
						endif
						if (xlop(n).lt.zero) then 
							pcello(n) = int(xlop(n)-1)
						else 
							pcello(n) = int(xlop(n))
						end if
					enddo
					!---------------------------------------

					!ON THE INNER SPHERE
					in1%theta(ith)%th = th
					in1%theta(ith)%phi(iphi)%ph = phi

					call interpolate_pdata(pcelli,xlip,pgli,pii,l)
					in1%theta(ith)%phi(iphi)%pres_comp(1) = pii
	!				in1%theta(ith)%phi(iphi)%pgrad = pgli
	!				call rotate_vector(in1%theta(ith)%phi(iphi)%pgrad,in1%theta(ith)%th,in1%theta(ith)%phi(iphi)%ph)

					call interpolate_udata(vcelli,xli,ib,ie,jb,je,kb,ke,uli,rdummy,rdummy,rdummy,0,m,l,onew)
					in1%theta(ith)%phi(iphi)%vel = uli
					call rotate_vector(in1%theta(ith)%phi(iphi)%vel,in1%theta(ith)%th,in1%theta(ith)%phi(iphi)%ph)

					!ON THE OUTTER SPHERE
					out1%theta(ith)%th = th
					out1%theta(ith)%phi(iphi)%ph = phi

					call interpolate_pdata(pcello,xlop,pglo,plo,l)
					out1%theta(ith)%phi(iphi)%pres_comp(1) = plo
	!				out1%theta(ith)%phi(iphi)%pgrad = pglo
	!				call rotate_vector(out1%theta(ith)%phi(iphi)%pgrad,out1%theta(ith)%th,out1%theta(ith)%phi(iphi)%ph)

					call interpolate_udata(vcello,xlo,ib,ie,jb,je,kb,ke,ulo,rdummy,rdummy,rdummy,0,m,l,onew)
					out1%theta(ith)%phi(iphi)%vel = ulo
					call rotate_vector(out1%theta(ith)%phi(iphi)%vel,out1%theta(ith)%th,out1%theta(ith)%phi(iphi)%ph)

	!				umagnormo  = dot_product(ulo(:),xx(:))
	!				do idim=1, ndim
	!					ulo_p(idim)  = ulo(idim)  - umagnormo  * xx(idim)
	!				enddo

					!ON THE SURFACE
					surf%theta(ith)%th = th
					surf%theta(ith)%phi(iphi)%ph = phi

					call interpolate_pdata(pcellb,xlp,pgl,pl,l)
					surf%theta(ith)%phi(iphi)%pres_comp(1) = pl
	!				surf%theta(ith)%phi(iphi)%pgrad = pgl
	!				call rotate_vector(surf%theta(ith)%phi(iphi)%pgrad,surf%theta(ith)%th,surf%theta(ith)%phi(iphi)%ph)

					call interpolate_udata(vcellb,xl,ib,ie,jb,je,kb,ke,ul,rdummy,rdummy,rdummy,0,m,l,onew) 
					surf%theta(ith)%phi(iphi)%vel = ul
					call rotate_vector(surf%theta(ith)%phi(iphi)%vel,surf%theta(ith)%th,surf%theta(ith)%phi(iphi)%ph)

	!				umagnorm   = dot_product(ul(:),xx(:))
	!				do idim=1, ndim
	!					ul_p(idim)   = ul(idim)   - umagnorm   * xx(idim)
	!				enddo
				enddo
			enddo

			do ith=0, nth
				th=ith*dth

				do iphi=0, nphi-1
					phi=iphi*dphi

					in1_vel  => in1%theta(ith)%phi(iphi)%vel
					surf_vel => surf%theta(ith)%phi(iphi)%vel
					out1_vel => out1%theta(ith)%phi(iphi)%vel

					in1_pres  => in1%theta(ith)%phi(iphi)%pres_comp(1)
					surf_pres => surf%theta(ith)%phi(iphi)%pres_comp(1)
					out1_pres => out1%theta(ith)%phi(iphi)%pres_comp(1)

	!				v_tmp = ufmean
	!				call rotate_vector(v_tmp,layer_c%th,layer_c%phi(iphi)%ph)
	!				v_tmp = layer_c%phi(iphi)%vel - v_tmp
	!				layer_c%phi(iphi)%vel = v_tmp

					vel_th(1,ith,:) = vel_th(1,ith,:) + in1_vel(:)
					pres_th(1,ith) = pres_th(1,ith) + in1_pres

					vel_th(2,ith,:) = vel_th(2,ith,:) + surf_vel(:)
					pres_th(2,ith) = pres_th(2,ith) + surf_pres

					vel_th(3,ith,:) = vel_th(3,ith,:) + out1_vel(:)
					pres_th(3,ith) = pres_th(3,ith) + out1_pres

					nullify(surf_vel, out1_vel, in1_vel)
					nullify(surf_pres, out1_pres, in1_pres)
				enddo
			enddo
		enddo

		do ith=0, nth
			vel_th(:,ith,:)     = vel_th(:,ith,:)     / nphi / (umeanslip/(1-maxvolfrac))
			pres_th(:,ith)      = pres_th(:,ith)      / nphi / (vis*(umeanslip/(1-maxvolfrac))**2)
		enddo

		if (I_AM_NODE_ZERO) then
			filename = trim(run_name)//"_"//trim(input_name)//"_surface_field.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")		

			do i=1, 3
				write (1,*) "zone"
				do ith=0, nth
					tmp1 = sqrt(dot_product(vel_th(i,ith,:),vel_th(i,ith,:)))
					write (1,"(1I6,6D15.7)") ith, ith*dth*180/pi, vel_th(i,ith,:), tmp1, pres_th(i,ith)
				enddo
			enddo
			close(1)
		endif
	end subroutine surface_field


	subroutine rotate_tensor(tensor,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: tensor(ndim,ndim)

		real(8) :: tmp(3,3)
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
				do dim3=1, ndim
					do dim4=1, ndim
						tmp(dim1,dim2) = tmp(dim1,dim2) + rotate(dim3,dim1) * rotate(dim4,dim2) * tensor(dim3,dim4)
					enddo
				enddo
			enddo
		enddo
		tensor(:,:) = tmp(:,:)
	end subroutine rotate_tensor

	subroutine rotate_vector(vec,th,phi)
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
	end subroutine rotate_vector

end module surface_module
