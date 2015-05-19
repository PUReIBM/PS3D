module random_num
#include "ibm.h"
	use precision  
	use global_data
	implicit none

#include "/home/mehr/SPRNG/include/sprng_f.h"

#if PARALLEL
#define USE_MPI 1
#endif

	integer :: streamnum1, streamnum2, streamnum3
	integer :: seed1, seed2, seed3
	integer :: nstreams, gtype
	SPRNG_POINTER :: stream1, stream2, stream3

contains
	subroutine init_rand
		implicit none
		character*100 filename

integer :: i

		gtype = 0

		streamnum1 = myid          !This stream is different on each proces
		streamnum2 = myid+nproc    !This stream is different on each proces
		streamnum3 = myid+nproc*2  !This stream is different on each proces
		nstreams  = nproc !extra stream is common to all processes
		seed1 = make_sprng_seed()
		seed2 = make_sprng_seed()
		seed3 = make_sprng_seed()
		stream1 = init_sprng(gtype,streamnum1,nstreams,seed1,SPRNG_DEFAULT)
		stream2 = init_sprng(gtype,streamnum2,nstreams,seed2,SPRNG_DEFAULT)
		stream3 = init_sprng(gtype,streamnum3,nstreams,seed3,SPRNG_DEFAULT)

		if (I_AM_NODE_ZERO) then
			filename = trim(run_name)//"_seed_used.dat"
			open (unit=1, file=trim(filename), status="replace", action="write")
			write (1,*) seed1, seed2, seed3
			close(1)
		endif
			

		write (*,*) myid, seed1, seed2, seed3



		do i=1, 1000000
			write (*,*) sprng(stream1), sprng(stream2), sprng(stream3)
		enddo
	end subroutine init_rand

	subroutine finalize_rand
		implicit none
		integer junk

		junk = free_sprng(stream1)
		junk = free_sprng(stream2)
		junk = free_sprng(stream3)
	end subroutine finalize_rand
end module random_nummber
