        subroutine ranjn(x,jf,jl,kf,kl,jdim,xm,cv,ndim)
c
c  routine to generate random vectors with a joint normal distribution
c  of mean xm and covariance matrix cv .
c
c  arguments:
c    x    - two-dimensional array in which random vectors are returned.
c            x(j,k) is the k th component of the random vector on the
c            j th trial.  the vector length is nk=kl+1-kf, with kf and
c            kl being the first and last values of k.  there are nj=jl+1-jf
c            trials, with jf and jl being the first and last trial numbers.
c            x has dimension x(jdim,kdim), where kdim .ge. kl .
c    xm   - vector of length n. xm(kk) is the mean of x(j,k), where
c            k=kf+kk-1.
c    cv   - two-dimensional array. cv(kk,ll) is the covariance of
c            x(j,k)and x(j,l), where k=kf+kk-1 and l=kf+ll-1.
c            (see further notes on cv below)
c    ndim - dimension of cv.
c
c  method:
c    the choleski decomposition of cv is performed so that
c    cv = l lt , where l is a lower triangular matrix and lt is
c    its transpose. then if y is a standardised normal random
c    vector, then  x = xm + l y  is a random vector with a joint
c    normal distribution of mean xm and covariance matrix cv.
c    subroutine rann2 is called to generate the standardised normal
c    random vectors.
c    cv is symmetric positive definite. only the lower triangle of
c    cv is used and is replaced by l. on the first call.
c    set cv(1,1) negative on first call to indicate that this cv has
c    already been decomposed.  on subsequent calls with the same cv,
c    the decomposition is then bypassed.
c
c  bugs:
c    a) if cv(1,1) < 0 on first call (invalid covariance matrix), ranjn
c        assumes that decomposition has already been performed and will
c        attempt to generate the random vectors with poor results.
c    b) if cv(1,1) = 0, decomposition will be done of every call (does
c        not lead to incorrect distributions).
c    c) there is no check that cv is symmetric - lower triangle only used.
c
c  other routines called - rann2
c
        dimension x(jdim,kl),xm(ndim),cv(ndim,ndim)
        save n
c
c  test for same covariance matrix
c
        if( cv(1,1) .lt. 0. ) then
         cv(1,1) = -cv(1,1)
         go to 130
         end if
c
c  check covariance matrix
c
        n=1+kl-kf
        do 50 j=1,n
        if( cv(j,j) .gt. 0. ) go to 50
        if( cv(j,j) .lt. 0. ) then
	   write(6,*) 'Stopped in ranjn.f'
	   write(6,*) 'Diagonal element less than zero'
	   goto 500
	endif
c
c  diagonal is zero : check that column is zero
c
        if( j .eq. n ) go to 50
        do 55 k=j+1,n
        if( cv(k,j) .ne. 0. ) then
	   write(6,*) 'Stopped in ranjn.f'
	   write(6,*) 'Diagonal element found to be zero'
	   write(6,*) 'Non zero element in column'
	   goto 500
	endif
55	continue
50      continue
c
c  perform choleski decomposition of cv using lower triangle.
c
        cv(1,1) = sqrt(cv(1,1))
        if( n .eq. 1 ) go to 130
c
        do 100 k=2,n
c
c  check that cv is a valid covariance matrix
c
          if( cv(1,1) .gt. 0. ) cv(k,1) = cv(k,1)/cv(1,1)
          sumllt = cv(k,1)**2
          if( k .eq. 2 ) go to 105
c
          do 110 j=2,k-1
            sum = 0.
c
            do 120 m=1,j-1
120        sum = sum+cv(j,m)*cv(k,m)
        if( cv(j,j) .gt. 0. ) cv(k,j) = (cv(k,j)-sum)/cv(j,j)
c
110      sumllt = sumllt+cv(k,j)**2
c
c  check that cv is positive definite
c
105     cvkk = cv(k,k)-sumllt
        if( cvkk .lt. 0. ) then
	   write(6,*) 'Stopped in ranjn.f'
	   write(6,*) 'Matrix non-positive definite'
	   goto 500
	endif
100     cv(k,k) = sqrt(cvkk)
c
130   continue
c
c  fill x with standardised normal random vectors
c
        call rann2(x,jf,jl,kf,kl,jdim)
c
c  rescale x
c
        do 200 j=jf,jl
c
c  multiply x by lower triangle of cv
c
          k = kl+1
          do 200 kk=n,1,-1
            k = k-1
c
c  sum row times x
c
            sum = 0.
            kd = kf-1
            do 210 ll=1,kk
              kd = kd+1
210        sum = sum+cv(kk,ll)*x(j,kd)
c
200      x(j,k) = xm(kk)+sum
c
        cv(1,1) = -cv(1,1)
        return
c
500   continue
        do 305 i=1,ndim
305     write(6,600) (cv(i,j),j=1,ndim)
        stop
c
600   format(1x,1p,20e10.3)
c
        end
