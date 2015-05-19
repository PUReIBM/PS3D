	subroutine rndg( x, n )
c
c  routine to generate discrete random numbers that approximate
c  a gaussian distribution
c
c  october 1990
c
c  x - vector in which random numbers are returned
c  n - number to be generated
c
c  other routines called:  rnu
c
c  method:  the random number returned, x, is a discrete random
c	variable, with equally-probable values g(i), i=1, nd.
c	once the discrete values g(i) have been formed, sampling
c	is trivial (see the 100 loop). 
c	the values of g(i) are formed as all possible sums of the
c	values of 4 other discrete variables, whose (positive) values,
c	given in data statements, are created by the program dggen.
c	the odd moments are zero by symmetry.  to machine accuracy
c	the variance is 1, and the flatness is 3.  with the
c	parameters set to n1=8, n2=9, n3=10, n4=12, the 6-th, 8-th
c	and 10-th moments are: 14.73, 98.30, and 811.6  (cf 15, 105,
c	945 for a gaussian).  for all x, the difference between the
c	cdf and the gaussian cdf (i.e. the error in the cdf) is less
c	than  0.005.
c
	parameter ( n1 = 8 , n2 = 9, n3 = 10 , n4 = 12 )
	dimension gc(n4,4), g(n1*n2*n3*n4), nc(4), amom(10)
	dimension x(n)
	double precision sum
	save nd,g
c
	data ifst /0/
	data ( gc(i,1), i = 1, 4 ) / 0.0791315883,  0.3056624234,
     1	                             0.6851569414,  1.8522603512  /
c
	data ( gc(i,2), i = 1, 4 ) / 0.1839743704,  0.4354168475,
     1	                             0.8178389072,  1.8993961811  /
c
	data ( gc(i,3), i = 1, 5 ) / 0.0856366232,  0.2853691578,
     1                  0.5440011621, 0.9214153290, 1.9406925440  /
c
	data ( gc(i,4), i = 1, 6 ) / 0.0801580697, 0.2553057969,
     1	 0.4600485563, 0.7161571383, 1.0778864622, 2.0104796886  /
c
c----  on first call, determine discrete values ----------------------
c
	if( ifst .eq. 0 ) then
	   ifst = 1
	   nd = n1*n2*n3*n4
	   nc(1) = n1
	   nc(2) = n2
	   nc(3) = n3
	   nc(4) = n4
c
c  complete component discrete gaussians
c
	   do 10 ig = 1, 4
	   nh = nc(ig) / 2
	   gc(nh+1,ig) = 0.
c
	   do 10 i = 1, nh
10	   gc( nc(ig)+1 - i , ig ) = - gc( i, ig )
c
c  form compound discrete gaussian
c
	   i = 0
	   do 50 j1 = 1, n1
	   g1 = gc( j1 , 1 )
c
	   do 50 j2 = 1, n2
	   g2 = g1 + gc( j2 , 2 )
c
	   do 50 j3 = 1, n3
	   g3 = g2 + gc( j3 , 3 )
c
	   do 50 j4 = 1, n4
	   i = i + 1
50	   g(i) = 0.5*( g3 + gc( j4 , 4 ) )
c
c  check g
c
	   icheck = 0
	   if( icheck .ne. 0 ) then
c
	      do 70 im = 1, 10
	      sum = 0.
c
	      do 80 i = 1, nd
80	      sum = sum + g(i) ** im
c
70	      amom(im) = sum / nd
c
	      write(0,*)' check of rndg '
	      write(0,600)amom
600	      format(1p,5e16.6)
	   endif
	endif
c
c-------------  end of set-up  ---------------------------------------
c
c  randomly sample from discrete values
c
	call rnu( x, n )
c
	do 100 j = 1, n
	i = 1 + x(j) * nd
100	x(j) = g(i)
c
	return
	end
