	subroutine vrand(i,x,s,n)
c
c	vrandr can be restarted from the calling pgm by setting
c	s = 0, which sets ifst=0.
c  modofied 8/13/88 by s.b. pope
c
c  routine to generate a vector of pseudo-random
c  numbers uniformly distributed on the exclusive
c  interval (0,1).
c  for use on 32 bit machines.
c
c  arguments:
c    i    integer seed ( i .ne. 0 ).
c    x    vector in which random numbers are returned.
c    n    dimension of vector x.
c    s    dummy argument - set to 1 in calling program for
c	    consistency with ap164 vrand.
c
c  method:
c    on the first call, a store of nstore random numbers is
c    generated.  random numbers from the store are then sampled
c    at random.  once sampled, the random number is immediately
c    replaced
c
c    for seed i (0 < i < 2**31-1), ix = mod(65539*i,2**31) is
c    a pseudo-random integer uniform on the exclusive interval
c    (0,2**31-1).  x = ix*0.46566127e-9 is then a pseudo-random
c    real number uniform on the exclusive interval (0.,1.).
c    i is replaced by ix on each call.
c    takes advantage of 32 bit machine's truncation of
c    integers >= 2**31.
c
c  notes:
c    in single precision (32 bit) 4.6566127e-10 <= x <= (1.0 - 5.96e-8).
c    if i = 0, seed locks at 0 and x(j) = 0. for all j.
c
c  other routines called - none.
c
	parameter( nstore = 97 )
	dimension store( nstore+1 )
	dimension x(n)
c
	save k,store
	data ifst/0/
c
c	restart if s=0.
c
	if(s.eq.0.)then
	  ifst=0
	  return
	endif
c
c  check for i = 0 - invalid seed
c
	if( i .eq. 0 ) go to 300
c
c  fill store on first call
c
	if( ifst .eq. 0 ) then
	   ifst = 1
c
	   do 10 j=1,nstore+1
	   i = i*65539
	   if( i .lt. 0 ) i=i+2147483647+1
	   xi = i
10	   store(j) = xi*0.46566127e-9
c
	   k = store( nstore+1) * nstore + 1
	endif

c
c  generate random numbers
c
	do 100 j=1,n
	 x(j) = store(k)
	 i = i*65539
	 if( i .lt. 0 ) i=i+2147483647+1
	 xi = i
	 store(k) = xi*0.46566127e-9
	 k = store(k) * nstore + 1
100    continue
	return
c
300   write(6,600)
	stop
c
600   format(/,'error in vrand - seed = 0.')
c
	end
