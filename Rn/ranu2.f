      subroutine ranu2(x,jf,jl,kf,kl,jdim)
c
c  ap164 compatible version of ranu2.
c  routine to generate vectors of independent random numbers
c  uniformly distributed in the exclusive interval (0,1).
c
c  arguments:
c    x - two dimensional array in which the random vectors are
c    returned.  x(j,k) is the kth component of the vector on
c    the jth trial.  the vector length is nk=kl+1-kf, with kf
c    and kl being the first and last values of k.  there are
c    nj=jl+1-jf trials, with jf and jl being the first and
c    last trial numbers.
c    x has dimensions x(jdim,kdim), where kdim .ge. kl.
c
c  notes:
c    for a one-dimensional array set jf=jl=1, or kf=kl=1.
c    a maximum of 8000 random numbers are generated at
c    each call to vrand; this limit is imposed by the
c    length of vector ran.
c    other routines called - vrand.
c
      common/rancom/ix/apran/ran(8000)
      dimension x(jdim,kl)
c
      nj = jl+1-jf
      nk = kl+1-kf
      nreq = nj*nk
      ngen = nreq
      if( ngen .gt. 8000 ) ngen=8000
      ncount = ngen
c
      call vrand(ix,ran,1,ngen)
c
      jk = 0
c
c  copy vector ran into 2d array x
c  generate more random numbers as needed
c
      do 100 k=kf,kl
      do 100 j=jf,jl
       jk = jk+1
       if( jk .le. 8000 ) go to 100
         jk = 1
         ngen = nreq-ncount
         if( ngen .gt. 8000 ) ngen=8000
         ncount = ncount+ngen
         call vrand(ix,ran,1,ngen)
c        write(10,800) ran(1),ran(ngen)
c800     format('ranu2: ran(1), ran(ngen) =',2i13)
100    x(j,k) = ran(jk)
c
      return
      end
