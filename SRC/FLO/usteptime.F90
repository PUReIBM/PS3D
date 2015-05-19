MODULE usteptime
#include "ibm.h"
	USE precision 
	USE constants 
	USE global_data
	use fftw3_interface
	implicit none 
Contains
	SUBROUTINE ustep
		IMPLICIT NONE 
		!-----------------------------------------------------------------------
		!	local variables
		REAL(prcn) :: b
		COMPLEX(prcn) ::  r, wtmp
		INTEGER :: i, j, k, idim

		!gauss_u = .true.
		do idim=1, ndim
#if !PARALLEL
			do k=1, local_no(3)
				do j=1, local_no(2)
					do i=1, local_no(1)
						if (w2(i,j,k)>small_number) then
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
						if(w2(i,j,k)>small_number) then
							if (idim == 1) then
								wtmp = wy(j)
							elseif (idim == 2) then
								wtmp = wz(k)
							elseif (idim == 3) then
								wtmp = wx(i)
							endif
#endif
							!-----------------------------------------------------------------------
							!	the left-hand-side
							b = one + half*dt*vis*w2(i,j,k)

							!-----------------------------------------------------------------------
							!	diffusion
							r = u(i,j,k,idim) * (one - half*dt*vis*w2(i,j,k))

							!-----------------------------------------------------------------------
							!	include convection terms
							!       sign change infront of dt done..
							!       to be consistent with Jamals convention
							!       remember nl term is negative..
							!       so negative sign, which should have been in front of dt,
							!       is now absorbed into the nl term
							r = r + dt * (coef(1,3)*nl(i,j,k,idim)+coef(1,4)*onl(i,j,k,idim)) + dt*ff(i,j,k,idim)*force_factor

							!-----------------------------------------------------------------------
							!	pressure terms
							r = r - dt * p(i,j,k)*wtmp

							u(i,j,k,idim) = r/b
						endif
					enddo
				enddo
			enddo
		enddo
		!new gauss_u = .false.
	END SUBROUTINE ustep
END MODULE usteptime
