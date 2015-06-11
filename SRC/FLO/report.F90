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

module report
#include "ibm.h"
	use global_data
	use dependent_functions, only : run_time_file_opener

contains

	subroutine report_actual_desired
		implicit none

		if (I_AM_NODE_ZERO) then
			WRITE(*,'(A25,3(2x,g17.8))')'USMEAN DES = ', usmean_des(:)
			WRITE(*,'(A25,3(2x,g17.8))')'USMEAN ACT = ', usmean_act(:)
			WRITE(*,'(A25,3(2x,g17.8))')'USMEAN MIX = ', usmean(:)
			WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN DES = ', ufmean_des(:)
			WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN ACTUAL = ', ufmean(:)
			!WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN CC = ', ufmean_cc(:)
			WRITE(*,'(A25,3(2x,g17.8))')'UMEAN ACTUAL = ', umean(:)
!			WRITE(*,'(A25,3(2x,g17.8))')'UMEAN TMP = ', umean_temp(:)
			WRITE(*,'(A40,3(2x,g12.5))') "MEANSLIPMOD AND REY(MEANSLIPMOD):", mixmeanslipmod, mixmeanslipmod*char_length /vis
		endif
	end subroutine report_actual_desired


	subroutine report_force
		implicit none

		real(prcn) :: norm_factor
		real(prcn) ::  norm_drag_chem, norm_drag_poly_spec(nphases), norm_drag_poly
		real(prcn) :: avg_force(ndim),avg_force_chem(ndim), avg_force_spec(nphases,ndim), avg_force_chem_spec(nphases,ndim), tmp_ferror_array(nerr_steps)
		integer :: pstart, pend, iphs, m, idim

		character*100 :: formfile, FILENAME1, FILENAME2, FILENAME3, FILENAME4, FILENAME5, FILENAME6, FILENAME7


		!-----------------------------------------------------------------------
		!      POST PROCESSING AND DIAGNOSTICS
		!-----------------------------------------------------------------------
		norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)*real(nbody,prcn)

		avg_force_spec = zero
		avg_force_chem_spec = zero
		force_chem = zero

		pstart = 1
		do iphs = 1, nphases
			pend = pstart + phase_array(iphs)%npart - 1
			do m=pstart, pend
				avg_force_spec(iphs,:) = avg_force_spec(iphs,:) + (pres(m,:)+visc(m,:)) - mpg(:)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0

				force_chem(m,:) =  (pres(m,:)+visc(m,:))

				avg_force_chem_spec(iphs,:) = avg_force_chem_spec(iphs,:) + force_chem(m,:)
			enddo
			avg_force_spec(iphs,:)      = avg_force_spec(iphs,:)/phase_array(iphs)%npart
			avg_force_chem_spec(iphs,:) = avg_force_chem_spec(iphs,:)/phase_array(iphs)%npart

			pstart = pend + 1
		end do

		avg_force = zero
		avg_force_chem = zero

		do iphs = 1, nphases
			avg_force(1:ndim) = avg_force(1:ndim) + phase_array(iphs)%volfrac*avg_force_spec(iphs,1:ndim)
			avg_force_chem(1:ndim) = avg_force_chem(1:ndim) + phase_array(iphs)%volfrac*avg_force_chem_spec(iphs,1:ndim)
		end do
    
		avg_force(1:ndim) = avg_force(1:ndim)/mean_volfrac
		avg_force_chem(1:ndim) = avg_force_chem(1:ndim)/mean_volfrac


!!$    do iphs = 1, nphases
!!$       if(I_AM_NODE_ZERO)WRITE(*,'(A,4(2x,g17.8))') 'FORCES: avg_force_spec : ', phase_array(iphs)%volfrac,avg_force_spec(iphs,1:ndim)
!!$    end do
!!$    if(I_AM_NODE_ZERO)WRITE(*,'(A,3(2x,g17.8))') 'FORCES: avg_force : ', avg_force(1:ndim)
		do iphs = 1, nphases
			norm_drag_spec(iphs) = DSQRT(avg_force_spec(iphs,1)**2.d0 + avg_force_spec(iphs,2)**2.d0 + avg_force_spec(iphs,3)**2.d0)
			norm_drag_chem_spec(iphs) = DSQRT(avg_force_chem_spec(iphs,1)**2.d0 + avg_force_chem_spec(iphs,2)**2.d0 + avg_force_chem_spec(iphs,3)**2.d0)
		end do
    
		norm_drag = DSQRT(avg_force(1)**2.d0 + avg_force(2)**2.d0 + avg_force(3)**2.d0)
		norm_drag_chem = DSQRT(avg_force_chem(1)**2.d0 + avg_force_chem(2)**2.d0 + avg_force_chem(3)**2.d0)
    

		mpg_avg = DSQRT(dot_product(mpg(1:ndim), mpg(1:ndim)))
		total_drag_mpg = mpg_avg*voldom/norm_factor
!		LHS_UF(:) = dufmeandt(:)*voldom -pres_total(:)/((one-maxvolfrac)) + visc_total(:)/((one-maxvolfrac))
!		MOD_LHS_UF = DSQRT(dot_product(LHS_UF(1:3), LHS_UF(1:3)))
!		MOD_LHS_UF = MOD_LHS_UF/norm_factor 
!		if(I_AM_NODE_ZERO)WRITE(*,'(A,5(2x,g17.8))') 'DRAG: PRES, VISC, TOTAL, MOD_LHS_UF/N, TOTAL_MPG', pres_drag,visc_drag, total_drag, MOD_LHS_UF, total_drag_mpg

		if(I_AM_NODE_ZERO)WRITE(*,'(A,5(2x,g17.8))') 'DRAG: PRES, VISC, TOTAL, TOTAL_MPG', pres_drag,visc_drag, total_drag, total_drag_mpg  
    
		if(TRIM(input_type).eq.'lubtest')then
			do iphs = 1, nphases
				norm_drag_spec(iphs) = norm_drag_spec(iphs)/(3.d0*pi*vis*(ucharmod)*phase_array(iphs)%dia)
				norm_drag_poly_spec(iphs) = norm_drag_chem_spec(iphs)
				norm_drag_chem_spec(iphs) = norm_drag_chem_spec(iphs)/(3.d0*pi*vis*(ucharmod)*phase_array(iphs)%dia)
				norm_drag_poly_spec(iphs) = norm_drag_poly_spec(iphs)!/(vis**2.d0)
			end do
			norm_drag  = norm_drag/(3.d0*pi*vis*(ucharmod)*char_length)
			norm_drag_poly  = norm_drag_chem
			norm_drag_chem  = norm_drag_chem/(3.d0*pi*vis*(ucharmod)*char_length)
			norm_drag_poly  = norm_drag_poly!/(vis**2.d0)
		else
			do iphs = 1, nphases
				norm_drag_spec(iphs) = norm_drag_spec(iphs)/(3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*phase_array(iphs)%dia)
				norm_drag_poly_spec(iphs) = norm_drag_chem_spec(iphs)
				norm_drag_chem_spec(iphs) = norm_drag_chem_spec(iphs)/(3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*phase_array(iphs)%dia)
				norm_drag_poly_spec(iphs) = norm_drag_poly_spec(iphs)!/(vis**2.d0)
			end do
			norm_drag  = norm_drag/(3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)
			norm_drag_poly  = norm_drag_chem
			norm_drag_chem  = norm_drag_chem/(3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)
			norm_drag_poly  = norm_drag_poly!/(vis**2.d0)
		end if

		do iphs = 1, nphases
			IF(norm_drag_spec(iphs).gt.ZERO)then
				phase_array(iphs)%ferror = ABS(norm_drag_spec(iphs) - phase_array(iphs)%fold)/norm_drag_spec(iphs)
			ELSE
				phase_array(iphs)%ferror = ONE
			END IF
		end do
		IF(norm_drag.gt.ZERO)THEN
			ferror = ABS(norm_drag - fold)/norm_drag
		ELSE
			ferror = one
		END IF
		fold = norm_drag
    
		do iphs = 1, nphases
			phase_array(iphs)%fold = norm_drag_spec(iphs)
		end do

!!$       IF(Re.EQ.Zero.and.ReT.gt.zero)Then
!!$          norm_drag1 = (norm_drag1)/(6.d0*pi*vis*SQRT(gran_temp)*radbdy(1)*dx)
!!$       ELSE
!!$          norm_drag1 = (norm_drag1)/(6.d0*pi*vis*ucharmod*radbdy(1)*dx)
!!$          norm_drag2 = norm_drag2/(6.d0*pi*vis*ucharmod*radbdy(nbody)*dx)
!!$          norm_drag  = norm_drag/(3.d0*pi*vis*ucharmod*dia_phys)
!!$       END IF

		IF(FROM_POST) RETURN
!		if(rks.eq.itrmax) then 
			!Rearrange the ferror_array array so that the last entry is flushed out
			tmp_ferror_array(1:nerr_steps) = ferror_array(1:nerr_steps)
			ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
			ferror_array(1) = ferror
			!PRINT*,'FERROR_A =', FERROR_ARRAY
			ferror_hist = SUM(ferror_array(1:nerr_steps))/nerr_steps
			do iphs = 1, nphases
				tmp_ferror_array(1:nerr_steps) = phase_array(iphs)%ferror_array(1:nerr_steps)
				phase_array(iphs)%ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
				phase_array(iphs)%ferror_array(1) = phase_array(iphs)%ferror
				!PRINT*,'FERROR_A =', FERROR_ARRAY
				phase_array(iphs)%ferror_hist = SUM(phase_array(iphs)%ferror_array(1:nerr_steps))/nerr_steps
			end do
!		end if
    
    
		IF(I_AM_NODE_ZERO.and.first_pass) THEN 
			FILENAME1 = TRIM(RUN_NAME)//'_norm_drag_chem'//'.dat'
			FILENAME2 = TRIM(RUN_NAME)//'_dragcoeffy'//'.dat'
			FILENAME3 = TRIM(RUN_NAME)//'_dragcoeffz'//'.dat'
			!!$ FILENAME2 = TRIM(RUN_NAME)//'_dragcoeffsum'//'.dat'
			FILENAME4 = TRIM(RUN_NAME)//'_norm_drag'//'.dat'
			FILENAME5 = TRIM(RUN_NAME)//'_force_part'//'.dat'
			FILENAME6 = TRIM(RUN_NAME)//'_drag_components'//'.dat'
			FILENAME7 = TRIM(RUN_NAME)//'_normdrag_poly'//'.dat'

			formfile='formatted'
			CALL  RUN_TIME_FILE_OPENER(unitnormdragchem,FILENAME1, formfile)
			CALL  RUN_TIME_FILE_OPENER(unitdragtavg,FILENAME2, formfile)
			CALL  RUN_TIME_FILE_OPENER(unitnormdrag,FILENAME4,formfile)
			CALL  RUN_TIME_FILE_OPENER(unitnormdragpoly,FILENAME7,formfile)
			CALL  RUN_TIME_FILE_OPENER(unitforce,FILENAME5,formfile)
			CALL  RUN_TIME_FILE_OPENER(unitdrag_comps,FILENAME6,formfile)
    	ENDIF

		if (I_AM_NODE_ZERO) then
			WRITE (*,'(A25,4(2x,g17.8))') 'NORM DRAGS, FERROR', norm_drag, ferror
			if(.not.((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'psd')))WRITE(*,'(A25,4(2x,g17.8))') 'NORM DRAGS PHASES:', norm_drag_spec(1:nphases)
!			IF(rks.eq.itrmax) then 
				!c_drag_st(1:nbody) = force(1:nbody,1)/(3.d0*pi*vis*dia_phys*ucharmod)
				WRITE (unitnormdragchem,'(500(2x, e20.12))') t/t_conv/(1-maxvolfrac), t/t_vis, t/t_diff, norm_drag_chem, norm_drag_chem_spec(1:nphases)

				!c_drag_st(1:nbody) = force(1:nbody,2)/(3.d0*pi*vis*dia_phys*ucharmod)
				!WRITE(unitdragtavg,'(15(2x,e20.12))')  t/t_conv/(1-maxvolfrac), t/t_vis, t/t_diff,(mean_drag_tavg(idim)/t,idim=1,ndim),(drag_tavg(1,idim)/t,idim=1,ndim), (drag_tavg(2,idim)/t,idim=1,ndim)
				!c_drag_st(1:nbody) = force(1:nbody,3)/(3.d0*pi*vis*dia_phys*ucharmod)
          
				WRITE(unitnormdrag,'(500(2x, e20.12))') t/t_conv/(1-maxvolfrac), t/t_vis, t/t_diff, norm_drag, norm_drag_spec(1:nphases), usmean_act(1:ndim), (ufmean(idim)/ucharmod,idim=1,ndim)
          
				WRITE(unitdrag_comps,'(20(2x, e20.12))') pres_drag, visc_drag, total_drag, total_drag_mpg, norm_drag, ABS(total_drag - total_drag_mpg)/(total_drag + SMALL_NUMBER)
          
				WRITE(unitnormdragpoly,'(500(2x, e20.12))') t/t_conv/(1-maxvolfrac), t/t_vis, t/t_diff, norm_drag_poly, norm_drag_poly_spec(1:nphases)
!			end IF
		end if
    
		! Checking for blown-up simulations
#if 0
		if(norm_drag.gt.1E+06) THEN 
			if (I_AM_NODE_ZERO) then 
				WRITE (*,'(A,2x,g17.8,2x,A)') 'NORM DRAG', norm_drag,' is greater than 1E+06'
				WRITE (*,'(A)') 'STOPPING THIS CASE AFTER WRITING THE BLOW UP INDICATOR FILE'
				OPEN (2000, file=TRIM(RUN_NAME)//'_BLOWUP.dat', form='formatted')

				close(2000, status="keep")
			end if
			PARALLEL_FINISH()
			STOP
		end if
#endif
	end subroutine report_force


end module report
