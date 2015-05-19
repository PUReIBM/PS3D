	subroutine rnu( x, n )
c
c  routine to generate pseudo-random numbers uniformly in the
c	exclusive interval [ 0 , 1 ].
c
c  november 1990
c
c   x(i) , i = 1, n  array in which random numbers are returned
c
c   the entries  rnuget  and  rnuput  can be used to get and set
c   the seeds i1 and i2.
c
c  method: see f. james, computer physics communications, vol.60
c	p.329, 1990.
c
	dimension x(n)
	save i1, i2
	data i1, i2 / 12345, 67890 /
c
	do 100 i = 1, n
c
	ik = i1 / 53668
	i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
	if( i1 .lt. 0 ) i1 = i1 + 2147483563
c
	ik = i2 / 52774
	i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
	if( i2 .lt. 0 ) i2 = i2 + 2147483399
c
	ix = i1 - i2
	if( ix .lt. 1 ) ix = ix + 2147483562
c
100	x(i) = float( ix ) * 4.656612e-10
c
	return
c
	entry rnuget( is1, is2 )
	is1 = i1
	is2 = i2
	return
c
	entry rnuput( is1, is2 )
	i1 = is1
	i2 = is2
	return
c
	end
