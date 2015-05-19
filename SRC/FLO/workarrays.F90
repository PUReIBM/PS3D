Module nlarrays 
	Use precision 
	Use constants 
	Use global_data
	implicit none

#if PARALLEL
	!REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: uatminus1, uatnxp2
#endif
end Module nlarrays

Module bcsetarrays
	Use precision 
	Use constants 
	Use global_data
	implicit none

	Real(prcn),  DIMENSION(:,:,:,:), ALLOCATABLE, target :: omega, fr, ppr, diffn
end Module bcsetarrays

Module nlmainarrays
	Use precision 
	Use constants 
	Use global_data
	implicit none

	real(prcn), DIMENSION(:,:,:,:), ALLOCATABLE,  Target ::  ubc, nlbc, onlbc
	real(prcn), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pbc

	real(prcn), Dimension(:,:,:,:), pointer ::  onlbcp, nlbcp, ubcp
	real(prcn), Dimension(:,:,:), pointer ::  pbcp
end Module nlmainarrays

module field_tmp_arrays
	use precision
	implicit none

	real(prcn), allocatable :: urtmp(:,:,:)
	complex(prcn), allocatable :: uftmp(:,:,:) !, uftmp2(:,:,:), uftmp3(:,:,:)
end module field_tmp_arrays
	



