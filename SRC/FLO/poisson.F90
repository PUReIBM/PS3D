module poisson
#include "ibm.h"
	USE precision 
	USE constants 
	USE global_data
	USE bcsetarrays
	use fftw3_interface

	IMPLICIT NONE 
	!///////////////////////////////////////////////////////////////////
	!	Calculate the pressure field using AB. The boundary conditions
	!	are such that only the flow field from the previous two
	!	timesteps (at the boundaries) is needed to solve the system.
	!	The pressure is estimated at (t+dt)
	!--------------------------------------------------------------------
contains 

	subroutine pressure
		implicit none 
		!-----------------------------------------------------------------------
		!	local variables
		complex(prcn) :: tmpc
		real(prcn) :: b
		complex(prcn) :: wtmp
		integer :: i, j, k, idim

#if !PARALLEL
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					if (w2(i,j,k)>small_number) then
						b = -w2(i,j,k)
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
		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					if (w2(i,j,k)>small_number) then
						b = -w2(i,j,k)
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

							tmpc = tmpc + wtmp * (coef(1,3)*nl(i,j,k,idim) + coef(1,4)*onl(i,j,k,idim) + ff(i,j,k,idim))
						enddo
						tmpc = tmpc / b
						p(i,j,k) = p(i,j,k) + prf * (tmpc-p(i,j,k))
					endif
				enddo
			enddo
		enddo
	end subroutine pressure

	subroutine divcorr
		implicit none 
		!-----------------------------------------------------------------------
		!	local variables
		real(prcn) ::  b
		complex(prcn) :: tmpc, wtmp, div
		integer :: i,j,k, idim

		!new gauss_phi= .true.
#if !PARALLEL
		do k=1, local_no(3)
			do j=1, local_no(2)
				do i=1, local_no(1)
					if (w2(i,j,k)>small_number) then
						b = -dt*w2(i,j,k)
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
		do k=1, local_no(2)
			do j=1, local_no(1)
				do i=1, local_no(3)
					if (w2(i,j,k)>small_number) then
						b = -dt*w2(i,j,k)
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
							!-----------------------------------------------------
							!	calculate divergence of intermediate velocity field
							!	form divergence in Fourier space on pressure grid
							tmpc = tmpc + wtmp * u(i,j,k,idim)
						enddo

						!---------------------------------
						!    velocity corection
						div = tmpc/b
						do idim=1, ndim
#if !PARALLEL
							if (idim == 1) then
								wtmp = wx(i)
							elseif (idim == 2) then
								wtmp = wy(j)
							elseif (idim == 3) then
								wtmp = wz(k)
							endif
#else
							if (idim == 1) then
								wtmp = wy(j)
							elseif (idim == 2) then
								wtmp = wz(k)
							elseif (idim == 3) then
								wtmp = wx(i)
							endif
#endif
							u(i,j,k,idim) = u(i,j,k,idim) - dt*wtmp*div
						enddo

						!---------------------------------
						!    pressure corection
						p(i,j,k) = p(i,j,k) + div + half*dt*vis*w2(i,j,k)*div
					endif
				enddo
			enddo
		enddo
		!new gauss_phi= .false.
	end subroutine divcorr
end module poisson
