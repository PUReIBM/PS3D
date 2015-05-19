	subroutine ranstd( x, n )
c
c  routine to standardize samples of random variables.
c
c  input:
c	x(i), i = 1, n	- samples of random variables
c	n		- number of random variables
c  output:
c	x(i), i = 1, n	- samples of random variables modified
c			  so that their sample mean and variance
c			  are zero and unity, respectively.
c	
	dimension x(n)
c
	if( n .le. 1 ) return
c
	sum1 = 0.
	sum2 = 0.
	do 100 i = 1, n
	sum1 = sum1 + x(i)
100	sum2 = sum2 + x(i)**2
c
	sum1 = sum1 / n
	sum2 = sum2 / n
c
	if( sum2 .gt. 0. ) then
	   sum2 = 1. / sqrt( sum2 - sum1**2 )
	else
	   write(6,*)' ranstd: n,sum2=',n,sum2,
     1			' cannot standardize variance'
	   sum2 = 1.
	endif
c
	sum1 = -sum1*sum2
	do 200 i = 1, n
200	x(i) = sum1 + sum2 * x(i)
c
	return
	end
