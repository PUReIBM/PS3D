module steptheflow
#include "ibm.h"
	!///////////////////////////////////////////////////////////////////////
	!	calculate velocity field at next time level
	!	pressure estimation and divergence correction done in place
	!	(for each wavenumber) to save storage
	!-----------------------------------------------------------------------
	USE precision 
	USE constants 
	Use poisson, Only : pressure, divcorr
	use usteptime 
	Use nlmainarrays, Only : pbcp, ubcp
	Use dependent_functions
	USE dem_mod

contains
	subroutine velstep(sflag,rks)
		use init_turb, only : check_divergence, calc_velreal
		use fftw3_interface
		use parallel

		implicit none

		integer :: j, k, i
		integer, intent(in) :: sflag, rks
		complex(prcn) :: usumloc, usum

		CALL pressure
		CALL ustep
		CALL divcorr
		if (debug_check) CALL check_divergence

		! COMMENTING THIS, BECAUSE GRAVITY IS CHANGED IN MPG. IF BOTH GRAVITY AND MPG ARE AVAILABLE, THEN UNCOMMENT THIS
		umean(:) = umean(:) + (-mpg(:) + frmean(:) - frame_accln(:))*dt !+ grav(:) 

		if(move_particles)then
			frame_vel(:) = frame_vel(:) + frame_accln(:)*dt
			frame_pos(:) = frame_pos(:) + frame_vel(:)*dt
		end if
    
		if (I_AM_NODE_ZERO) then
			WRITE(*,'(A25,3(2x,g17.8))')'FR MEAN    @ N  = ', (FRMEAN(j), j = 1, ndim)
			WRITE(*,'(A25,3(2x,g17.8))')'MPG        @ N  = ', (MPG(j), j = 1, ndim)
			if (move_particles) then
				WRITE(*,'(A25,3(2x,g17.8))')'FRAME ACCLN@ N  = ', (frame_accln(j), j = 1, ndim)
				WRITE(*,'(A25,3(2x,g17.8))')'FRAME VEL  @ N  = ', (frame_vel(j), j = 1, ndim)
				WRITE(*,'(A25,3(2x,g17.8))')'FRAME POS  @ N  = ', (frame_pos(j), j = 1, ndim)
			endif
			WRITE(*,'(A25,3(2x,g17.8))')'UMEAN      @ N+1= ', (UMEAN(j), j = 1, ndim)
		endif
	end subroutine velstep
end module steptheflow
