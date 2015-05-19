      subroutine ranchi(x,jf,jl,kf,kl,jdim)
c
c  routine to generate vectors of independent random numbers with
c  a user-specified distribution. - Chi squared.
c
c  **  set up for pdf  fdens(x) = sqrt(chi squared) function
c
c  arguments:
c       x - two-dimensional array in which the random numbers are
c       returned.  x(j,k) is the k th component of the random vector
c       on the j th trial.  the vector length is nk=kl+1-kf, with kf
c       and kl being the first and last values of k.  there nj=jl+1-jf
c       trials , with jf and jl being the first and last trial numbers.
c       x has dimensions x(jdim,kdim), where kdim .ge. kl .
c irtn is a return code - set irtn=1 in case of convergence
c failure in table setup.
c
c  notes:
c       user must supply functions - see below.
c       for one-dimensional arrays, jf=jl=1.
c       other routines called - ranu2
c
c  method:
c       the random variable x has the specified distribution function
c       fdist(x) and the density function fdens(x).  if y is a uniformly
c       -distributed random variable (0. lt. y .lt. 1.), then x
c       satisfying the implicit equation  y=fdist(x)  is the required
c       random number.  on the first call, a table of x values is set
c       up corresponding to equally spaced y values.  then the
c       subroutine ranu2 is called to obtain random values of y, and
c       the corresponding x value is obtained from the table by
c       linear interpolation.  if the inverse of the distribution
c       function, finv, is known, then the table can be set up by
c         x = finv( fdist(x) ) = finv(y) .    in this case the user
c       must supply the fortran statement function finv(y).
c       if the inverse is not known then the distribution function
c       fdist(x) and the pdf fdens(x) must be supplied as statement
c       functions.  an iterative procedure is then used to solve the
c       implicit equation  y = fdist(x) .
c
c  user specifications:
c       ntab - number of table entries.
c       mode=1 - inverse function finv is supplied,
c           =2 - fdist and fdens are supplied.
c       xmin,xmax - minimum and maximum values of x. (even if x is
c               unbounded, finite values must be specified)
c       finv(y) - inverse of fdist(x) - mode=1 only.
c       fdist(x) - distribution function. it must be strictly increasing
c               and be between zero and one - mode=2 only.
c       fdens(x) - pdf (derivative of fdist(x)).  it must be strictly
c               positive everywhere - mode=2 only.
c	ndeg - number of degrees of freedom for the chi-squared.
c
c  comments:
c       this routine is designed to be computationally efficient,
c       with moderate accuracy, for the generation of a large number
c       of random numbers.  it is not extremely accurate and is
c       inefficient for generating only a few numbers.
c
      dimension tab(201),dtab(201),x(jdim,kl)
      save tab,dtab,ntabm,init
      data ntab/201/
      data xmin,xmax,mode/ 0.0, 3.29 , 2 /
      data ave, sd / 2.0, 2.0 /
      data init/0/
c
      if(init .eq. 0 ) go to 500
c
c  fill array with equally spaced y values.
c
c50    call ranu2(x,jf,jl,kf,kl,jdim)
50    dy=1./float(jl-jf+1)
      do 98 k=kf,kl
	x(jf,k) = dy
	x(jl,k) = 1.
        y=dy
        do 99 j=jf+1,jl-1
	  y = y + dy 
	  x(j,k)=y
  99	continue
  98   continue
c
c  linearly interpolate from table
c
      do 100 j=jf,jl
      do 100 k=kf,kl
      y=x(j,k)*ntabm+1
      i=y
      y=y-i
      x(j,k)=tab(i)+y*dtab(i)
100   continue
c
      return
c
c  set up table
c
500   init=1
      ntabm=ntab-1
      dy=1./float(ntabm)
      tab(1)=xmin
      tab(ntab)=xmax
c
c  set up table from inverse of fdist
c
      if( mode .eq. 2 ) go to 600
      y=0.
      do 510 i=2,ntabm
      y=y+dy
510   tab(i)=finv(y)
      go to 700
c
c  set up table by solving y=fdist(x).  combination of bisection and
c  newton,s method.
c
600   y=0.
      xx=xmin+dy*(xmax-xmin)
      do 610 i=2,ntabm
      y=y+dy
c set max. and min. xx and initial guess
      xxmin=tab(i-1)
      xxmax=xmax
c  iterate
      do 620 iter=1,50
      yy=fdist(xx,ave,sd)
      err=yy-y
      if( abs(err) .lt. 1.e-6 ) go to 610
      if( err .gt. 0. ) xxmax=xx
      if( err .lt. 0. ) xxmin=xx
c  newton,s method
      xx=xx+err/amax1(fdens(xx,ave,sd),1.e-10)
c if new value exceeds known bounds, use bisection
      if( xx .le. xxmin .or. xx .ge. xxmax ) xx=0.5*(xxmin+xxmax)
620   continue
c
c  convergence failure - return to calling program
c
      irtn = 1
      return
c
c  activate the following statements for debugging in case of
c  convergence failure
c
c     write(6,900) y,err,xx,xxmin,xxmax
c900  format(1x,'convergence failure in ranal2',1p5e13.3)
c     stop
c
c  convergence - enter value in table
c
610   tab(i)=xx
c
c  set dtab
c
700   do 710 i=1,ntabm
710   dtab(i)=tab(i+1)-tab(i)
c
c  end of table set up
c
      go to 50
c
      end
c
c now for the functions
c	ndeg = number of degrees of freedom.
c	ave = ndeg = a/b
c	sd = ndeg*2 = a/b**2
c	a = n/2
c	b = 1/2
c
      function fdist( xx, ave, sd )
      a = ave*0.5
      z = xx**2
      fdist = gammp(a,z*0.5)
      return
      end
c
      function fdens( xx, ave, sd )
      z = xx**2
      a = ave*0.5
      glnb2 = exp(gammln(a))
      fdens = 0.5*xx*z**(a-1.0) * ( exp(-z*0.5) ) / ((2.0**a) * glnb2)
      return
      end
c
      function finv( xx, ave, sd )
      finv = 0.0
      return
      end
c
c
c
      function gammln(xx)
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end
c
      function gammp(a,x)
c
c	Returns the incomplete gamma function P(a,x) - see Num Recipes.
c
      if(x.lt.0..or.a.le.0.)pause
      if(x.lt.a+1.)then
	call gser(gamser,a,x,gln)
	gammp = gamser
      else
	call gcf(gammcf,a,x,gln)
	gammp = 1.-gammcf
      endif
      return
      end
c
c
	subroutine gser(gamser,a,x,gln)
c
c	returns the incomplete gamma fn P(a,x) evaluated by its series
c	representation as GAMSER. Also returns ln Gamma(a) as gln.
c	
	parameter (itmax=100,eps=3.e-7)
	gln=gammln(a)
	if(x.le.0.)then
	  if(x.lt.0.)pause
	  gamser=0.
	  return
	endif
	ap=a
	sum=1./a
	del=sum
	do 11 n=1,itmax
	  ap=ap+1.
	  del=del*x/ap
	  sum = sum + del
	  if(abs(del).lt.abs(sum)*eps)goto 1
  11	continue
	pause 'A too large, ITMAX too small'
  1	gamser = sum*exp(-x+a*log(x)-gln)
	return
	end
c
c
c
	subroutine gcf(gammcf,a,x,gln)
c
c	returns the incomplete gamma function Q(a,x) evaluated by its
c	continued fraction representation as GAMMCF. Also returns
c	Gamma(a) as GLN.
	parameter (itmax=100,eps=3.e-7)
	gln = gammln(a)
	gold=0.
	a0=1.
	a1=x
	b0=0.
	b1=1.
	fac=1.
	do 11 n=1,itmax
	  an=float(n)
	  ana= an-a
	  a0 = (a1+a0*ana)*fac
	  b0 = (b1+b0*ana)*fac
	  anf  = an*fac
	  a1 = x*a0+anf*a1
	  b1 = x*b0+anf*b1
	  if(a1.ne.0.)then
	    fac=1./a1
	    g = b1*fac
	    if(abs((g-gold)/g).lt.eps)goto 1
	    gold=g
	  endif
  11	continue
	pause 'A too large, ITMAX too small'
 1	gammcf=exp(-x+a*alog(x)-gln)*g
	return
	end
      function betacf(a,b,x)
      parameter (itmax=100,eps=6.e-6)
      am=1.
      bm=1.
      az=1.
      qab=a+b
      qap=a+1.
      qam=a-1.
      bz=1.-qab*x/qap
      do 11 m=1,itmax
        em=m
        tem=em+em
        d=em*(b-m)*x/((qam+tem)*(a+tem))
        ap=az+d*am
        bp=bz+d*bm
        d=-(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
        app=ap+d*az
        bpp=bp+d*bz
        aold=az
        am=ap/bpp
        bm=bp/bpp
        az=app/bpp
        bz=1.
        if(abs(az-aold).lt.eps*abs(az))go to 1
11    continue
      pause 'a or b too big, or itmax too small'
1     betacf=az
      return
      end
