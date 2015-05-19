Module useful_functions

#include "../FLO/ibm.h"

  USE precision  
  USE constants  
  USE nr
  USE nrutil
  !Author: Rahul Garg
  Implicit None 
  INTERFACE mis_avg_data
     MODULE PROCEDURE MIS_AVG_REAL_ZEROD_DATA
     MODULE PROCEDURE MIS_AVG_REAL_ONED_DATA
     !     MODULE PROCEDURE MIS_AVG_INTEGER_ONED_DATA
  END INTERFACE
Contains 
  SUBROUTINE MIS_AVG_REAL_ZEROD_DATA(data, avg_data, err_bar, conf1, mis_true_in, CHECK_FOR_OUTLIER)
    Implicit none 
    REAL(prcn), DIMENSION(:), INTENT(INOUT) :: data 
    REAL(prcn) , INTENT(OUT) :: avg_data, err_bar
    REAL(prcn), intent(in) :: conf1
    INTEGER, DIMENSION(:), INTENT(IN) :: mis_true_in
    INTEGER, PARAMETER :: SP = KIND(1.0)
    LOGICAL, OPTIONAL  :: CHECK_FOR_OUTLIER
    INTEGER, DIMENSION(size(mis_true_in,1)) :: mis_true
    LOGICAL :: data_delete
    INTEGER :: ndata, idata, data_count, delete_count, data_pts_orig
    REAL(prcn) :: dvar, sdev, ctemp, tmp, ds, xminmubysig 
    REAL(SP) :: prob_xminmubysig, ntimesprob

    REAL(PRCN) :: confint 
    delete_count = 0
    data_pts_orig = SUM(mis_true_in)
    if(.not.present(check_for_outlier)) check_for_outlier = .false.

    !WRITE(*,*) 'CHECK FOR OUTLIER ?', CHECK_FOR_OUTLIER
    mis_true = mis_true_in
2000 continue 
    ndata = SIZE(mis_true,1) 

    tmp = zero 
    data_count = 0
    do idata = 1, ndata 
       if(mis_true(idata).eq.1) then 
          tmp = tmp + data(idata)
          data_count = data_count + 1
       end if
    end do

    avg_data = tmp/real(data_count)
    dvar=0.d0
    DO idata=1,ndata
       if(mis_true(idata).eq.1) then 
          ds=data(idata)-avg_data
          dvar=dvar+ds*ds
       end if

    ENDDO


    IF(data_count.gt.1) then 
       CALL GET_CONFIN(data_count, confint)

       dvar=dvar/DBLE(float(data_count-1))

       sdev=DSQRT(dvar)

       ctemp = sdev/SQRT(float(data_count))
       err_bar = confint*ctemp
    ELSE
       err_bar = zero 
    end IF

    data_delete = .false.

    if(check_for_outlier) then 
       IF(data_count.ge.5) then 
          do idata = 1, ndata 
             if(mis_true(idata).eq.1) then 
                xminmubysig = ABS((data(idata)-avg_data)/sdev)
                prob_xminmubysig = erfc( real(xminmubysig,sp))
                ntimesprob = real(data_count,prcn)*real(prob_xminmubysig,prcn)
                if(ntimesprob.lt.0.5d0) then 
                   WRITE(*,'(A,3(2x,g17.8))')'MU, SIG, X-MU/SIG = ', avg_data, sdev, xminmubysig
                   WRITE(*,'(A,g17.8)')'rejecting', data(idata)
                   WRITE(*,'(A,100(/,2x,g17.8))')'OVERALL DATA ', data(1:data_pts_orig)
                   !WRITE(*,'(A,100(/,2x,g17.8))')'OVERALL DATA ', data(1:data_pts_orig)
                   mis_true(idata) = 0
                   data_delete = .true. 
                   delete_count = delete_count + 1
                   WRITE(*,'(A,2x,i4)')'POINTS DELETED SO FAR = ', delete_count
                   goto 1000
                end if

             end if
          end do
       end IF
    end if
    
1000 continue 
    if(data_delete) goto 2000
  end SUBROUTINE MIS_AVG_REAL_ZEROD_DATA

  SUBROUTINE MIS_AVG_REAL_ONED_DATA(data, avg_data, err_bar, conf1, mis_true)!, check_for_outlier)
    Implicit none 
    REAL(prcn), DIMENSION(:,:), INTENT(INOUT) :: data 
    REAL(prcn) , DIMENSION(:), INTENT(OUT) :: avg_data, err_bar
    REAL(prcn), intent(in) :: conf1
    INTEGER, DIMENSION(:) :: mis_true
    INTEGER :: ndata, idata, datalen, idatalen
    !LOGICAL, OPTIONAL  :: CHECK_FOR_OUTLIER
    LOGICAL :: check_for_outlier2

    ndata = SIZE(mis_true,1)

    datalen = SIZE(data,1)

    do idatalen = 1, datalen

       CALL MIS_AVG_REAL_ZEROD_DATA(data(idatalen, 1:ndata), avg_data(idatalen), err_bar(idatalen), conf1, mis_true(1:ndata), check_for_outlier = .false.)
    end do

  end SUBROUTINE MIS_AVG_REAL_ONED_DATA

  SUBROUTINE initialize_lookup_table(tunit)
    USE precision  
    USE constants  
    USE general_funcs
    USE post_global_data

    IMPLICIT NONE

    CHARACTER*8 :: ciperc(percmax)
    INTEGER, Intent(in) :: tunit
    INTEGER :: iperc,isamples

    READ(tunit,*)(ciperc(iperc),iperc=1,percmax)

    do isamples=1,nmismax
       cis(isamples)%confper(1:percmax) = ciperc(1:percmax)
    end do
    do isamples=1,nmismax
       READ(tunit,*)(cis(isamples)%confi(iperc),iperc=1,percmax)
    end do

  END SUBROUTINE initialize_lookup_table

  RECURSIVE function lookup_entry(nsamples,entry) result(confis)
    USE precision  
    USE constants  
    USE post_global_data
    IMPLICIT NONE
    INTEGER, Intent(in) :: nsamples
    CHARACTER*8, Intent(inout) :: entry
    REAL(prcn) :: confis
    Integer :: iperc

    do iperc = 1,percmax
       if(entry.eq.(cis(nsamples)%confper(iperc)))then
          confis = cis(nsamples)%confi(iperc)
          exit
       end if
    end do

    if(iperc.gt.percmax)then
       PRINT*,'ENTRY ', entry,' NOT FOUND. SETTING 95% CONFIDENCE INTERVAL'
       entry = TRIM('95%')
       confis = lookup_entry(nsamples,entry)
       RETURN
    end if
    PRINT*,'CONFIDENCE INTERVAL FOR ', TRIM(entry),' CONFIDENCE AND ',nsamples, ' DEGREES OF FREEDOM IS ', confis 
    RETURN
  END function lookup_entry

  SUBROUTINE GET_CONFIN(ndatain, confin)
    USE post_global_data, Only : nmismax, per_conf 
    USE global_data, Only : minunitno, maxunitno 
    USE general_funcs
    implicit none 
    INTEGER, INTENT(IN) :: ndatain 
    REAL(prcn), INTENT(OUT) ::  confin
    CHARACTER*80 :: table_name
    INTEGER :: tunit, ndata 
    LOGICAL :: filexist,isopen
    ndata = ndatain 
    !  WRITE(*,*)'IN GET CONFIN, NDATA = ', ndata
    table_name = "ci_lookup.dat"
    INQUIRE(FILE=table_name,EXIST=filexist,OPENED=isopen)
    if(filexist)then
       tunit = getnewunit(minunitno,maxunitno)
       OPEN(unit=tunit,FILE=table_name,form='formatted')
       CALL initialize_lookup_table(tunit)
       if(ndata.gt.nmismax)ndata=nmismax
       if(ndata.gt.1)then
          confin = lookup_entry(ndata-1,per_conf)

          !WRITE(*,'(2(A,2x,g17.8))') 'DATA PTS = ', real(NDATA),' CONFIN = ', CONFIN
       else
          !WRITE(*,'(A,2x,g17.8)') 'ONLY ONE MIS: CONFIN = ', CONFIN
          confin = zero
       end if
       close(tunit,status='keep')

    else

       confin = 1.812
    end if
  end SUBROUTINE GET_CONFIN

#if 0
  SUBROUTINE jacobi(a,n,np,d,v,nrot)
    IMPLICIT NONE
    Integer,Intent(in):: n,np
    Integer :: NMAX,nrot
    double precision a(np,np), d(np), v(np,np)
    PARAMETER (NMAX=500)
    !      Computes all eigenvalues and eigenvectors of a real symmetric N ×N
    !matrix a. On output,
    !      elements of a above the diagonal are destroyed. d is a vector of length
    !N that returns the
    !      eigenvalues of a. v is an N ×N matrix whose columns contain, on output,
    !the normalized
    !      eigenvectors of a. nrot returns the number of Jacobi rotations that were
    !required.
    INTEGER i,ip,iq,j
    double precision c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
    do ip=1,n
       do iq=1,n
          v(ip,iq) = 0.0
       enddo
       v(ip,ip) = 1.0
    enddo
    !      do ip=1,n
    !         print*,(a(ip,iq),iq=1,2)
    !      enddo
    do ip=1,n
       b(ip)= a(ip,ip)
       d(ip) = b(ip)
       !         print*,ip,d(ip),a(ip,ip)
       z(ip) = 0.0
    enddo
    nrot=0
    do i=1,500
       !sum = 0.0
       do ip=1,n-1
          do iq=ip+1,n
             sm=sm+abs(a(ip,iq))
          enddo
       enddo
       if (sm.eq.0.0) then
          print*,'return'
          return
       endif
       if (i.lt.4) then
          tresh=0.2*sm/n**2
       else
          tresh =0.0
       endif
       do ip=1,n-1
          do iq=ip+1,n
             g=100.0*abs(a(ip,iq))
             if ((i.gt.4).and.(abs(d(ip))&
                  &+g.eq.abs(d(ip))).and.(abs(d(iq))&
                  &+g.eq.abs(d(iq)))) then
                a(ip,iq) = 0.0
             else if(abs(a(ip,iq)).gt.tresh)then
                h=d(iq)-d(ip)
                if (abs(h)+g.eq.abs(h))then
                   t=a(ip,iq)/h
                else
                   theta=0.5*h/a(ip,iq)
                   t=1.0/(abs(theta)+sqrt(1.0+theta**2))
                   if (theta.lt.0.0) t=-t
                endif
                c=1.0/sqrt(1.0+t**2)
                s=t*c
                tau = s/(1.0+c)
                h=t*a(ip,iq)
                !                  print*,d(1),d(2),h
                z(ip) = z(ip)-h
                z(iq) = z(iq)+h
                d(ip) = d(ip)-h
                d(iq) = d(iq)+h
                a(ip,iq) = 0.0
                do j=1,ip-1
                   g=a(j,ip)
                   h=a(j,iq)
                   a(j,ip) = g-s*(h+g*tau)
                   a(j,iq) = h+s*(g-h*tau)
                enddo
                do j=ip+1,iq-1
                   g=a(ip,j)
                   h=a(j,iq)
                   a(ip,j) = g-s*(h+g*tau)
                   a(j,iq) = h+s*(g-h*tau)
                enddo
                do j=iq+1,n
                   g=a(ip,j)
                   h=a(iq,j)
                   a(ip,j) = g-s*(h+g*tau)
                   a(iq,j) = h+s*(g-h*tau)
                enddo
                do j=1,n
                   g=v(j,ip)
                   h=v(j,iq)
                   v(j,ip) = g-s*(h+g*tau)
                   v(j,iq) = h+s*(g-h*tau)
                enddo
                nrot=nrot+1
                !                  print*,i,nrot,h
                !                  print*,d(1),d(2)
                !                  print*,a(1,1),a(1,2)
                !                  print*,a(2,1),a(2,2)
             endif
          enddo
       enddo
       do ip=1,n
          b(ip) = b(ip)+z(ip)
          d(ip) = b(ip)
          z(ip) = 0.0
       enddo
    enddo
    !      pause 'too many iterations in jacobi'
    return
  END SUBROUTINE jacobi

  subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
!!$c
    integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
    double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
    double precision c,f,g,h,s,x,y,z,scale,anorm
    logical matu,matv
!!$c
!!$c       This subroutine is a translation of the ALGOL procedure
!!$c       SVD Num. Math. 14,403-420(1970) by Golub and Reinsch.
!!$c       Handbook for Auto. Comp., Vol II-Linear Algebra, 134-151(1971).
!!$c
!!$c       This subroutine determines the singular value decomposition
!!$c        T
!!$c       A = USV of a real M by N rectangular matrix. Householder
!!$c       bidiagonalization and a variant of the QR algorithm are used.
!!$c
!!$c       On input
!!$c
!!$c         NM must be set to the row dimension of two-dimensional array
!!$c            parameters as declared in the calling program dimension
!!$c            statement. Note that nm must be atleast as large of the
!!$c            maximum of M and N.
!!$c
!!$c         M is the number of rows of A ( and U).
!!$c
!!$c         N is the number of columns of A ( and U) and the order of V.
!!$c
!!$c         A contains the rectangular input matrix to be decomposed.
!!$c
!!$c         MATU should be set to .TRUE. if the U matrix in the decomposition
!!$c              is desired., and to .FALSE. otherwise.
!!$c
!!$c
!!$c         MATV should be set to .TRUE. if the V matrix in the decomposition
!!$c              is desired., and to .FALSE. otherwise.
!!$c
!!$c       On output.
!!$c
!!$c       A is unaltered ( unless overwritten by U or V).
!!$c
!!$c       W contains the N(non-negative) singular values of A ( the diagonal
!!$c         elements of S ). They are unordered. If an error exit is made then
!!$c         the singular values should be correct for indices ierr+1, ierr+2,..N
!!$c
!!$c       U contains the matrix U (orthogonal column vectors) of the decompo-
!!$c         sition if MATU has been set to .TRUE. Otherwise U is used as a
!!$c         temporary array.U may coincide with A. If an error exit is made,
!!$c         the columns of U corresponding to indices of correct singular values
!!$c         should be correct.
!!$c
!!$c       V contains the matrix V (orthogonal column vectors) of the decompo-
!!$c         sition if MATV has been set to .TRUE. Otherwise V is used as a
!!$c         temporary array.V may coincide with A. If an error exit is made,
!!$c         the columns of V corresponding to indices of correct singular values
!!$c         should be correct.
!!$c
!!$c       IERR is set
!!$c        ZERO for normal return,
!!$c        K    if the K-th singular value has not been determined after
!!$c             30 iterations.
!!$c
!!$c        RV1 is a temporary storage array.
!!$c
!!$c       Questions and comments should be directed to B. S. Garbow, Applied
!!$c       Mathematics Divisions, Argonne National Laboratory.
!!$c
!!$c       Modified to eliminate MACHEP.
!!$c
    ierr=0

    do 100 i = 1, m

       do 100 j = 1, n
          u(i,j) = a(i,j)
100       continue
!!$c       .... Householder reduction to bidiagonal form ......
          g = 0.0
          scale = 0.0
          anorm = 0.0

          do 300 i = 1, n
             l = i + 1
             rv1(i) = scale * g
             g = 0.0
             s = 0.0
             scale = 0.0
             if(i .gt. m) goto 210

             do 120 k = i, m
120             scale = scale + abs(u(k,i))

                if(scale.eq.0.0) goto 210

                do 130 k = i, m
                   u(k,i) = u(k,i)/scale
                   s = s + u(k,i)**2
130                continue

                   f = u(i,i)
                   g = -sign(sqrt(s),f)
                   h = f * g - s
                   u(i,i) = f - g
                   if( i .eq. n) goto 190

                   do 150 j = l, n
                      s = 0.0

                      do 140 k = i, m
140                      s = s + u(k,i) * u(k,j)

                         f = s / h

                         do 150 k =  i, m
                            u(k,j) = u(k,j) + f * u(k,i)
150                         continue

190                         do 200 k = i, m
200                            u(k,i) = scale * u(k,i)

210                            w(i) = scale * g
                               g = 0.0
                               s = 0.0
                               scale = 0.0
                               if(i.gt.m .or. i .eq. n) goto 290

                               do 220 k = l, n
220                               scale = scale + abs(u(i,k))

                                  if(scale .eq. 0.0) goto 290

                                  do 230 k = l, n
                                     u(i,k) = u(i,k) / scale
                                     s = s + u(i,k)**2
230                                  continue

                                     f = u(i,l)
                                     g = -sign(sqrt(s),f)
                                     h = f * g - s
                                     u(i,l) = f - g

                                     do 240 k = l, n
240                                     rv1(k) = u(i,k)/ h

                                        if(i .eq. m) goto 270

                                        do 260 j = l, m
                                           s = 0.0

                                           do 250 k = l, n
250                                           s = s + u(j,k) * u(i,k)

                                              do 260 k = l, n
                                                 u(j,k) = u(j,k) + s * rv1(k)
260                                              continue

270                                              do 280 k = l, n
280                                                 u(i,k) = scale * u(i,k)

290                                                 anorm = dmax1(anorm,(dabs(w(i))+dabs(rv1(i))))
300                                                 continue
!!$c       ....accumulation of right-hand transformations .......
                                                    if(.not.matv) goto 410
!!$c       ....for i = n step -1 until 1 do -- ........
                                                    do 400 ii = 1, n
                                                       i = n + 1 - ii
                                                       if(i .eq. n) goto 390
                                                       if(g .eq. 0.0) goto 360

                                                       do 320 j = l, n
!!$c       ....double division avoids possible underflow .....
320                                                       v(j,i) = ( u(i,j) / u(i,l) ) / g

                                                          do 350 j = l, n
                                                             s = 0.0

                                                             do 340 k = l, n
340                                                             s = s + u(i,k) * v(k,j)

                                                                do 350 k = l, n
                                                                   v(k,j) = v(k,j) + s * v(k,i)
350                                                                continue
360                                                                do 380 j = l, n
                                                                      v(i,j) = 0.0
                                                                      v(j,i) = 0.0
380                                                                   continue

390                                                                   v(i,i) = 1.0
                                                                      g = rv1(i)
                                                                      l = i
400                                                                   continue
!!$c ....accumulation of left-hand transformations .....
410                                                                   if(.not.matu) goto 510
!!$c .... for i=min(m,n) step -1 until 1 do -- ......
                                                                      mn = n
                                                                      if ( m .lt. n ) mn = m

                                                                      do 500 ii = 1, mn
                                                                         i = mn + 1 - ii
                                                                         l = i + 1
                                                                         g = w(i)
                                                                         if(i .eq. n) goto 430

                                                                         do 420 j = l, n
420                                                                         u(i,j) = 0.0

430                                                                         if(g.eq.0.0) goto 475
                                                                            if(i.eq.mn) goto 460

                                                                            do 450 j = l, n
                                                                               s = 0.0

                                                                               do 440 k = l, m
440                                                                               s = s + u(k,i) * u(k,j)
!!$c   ....double division avoids possible underflow ....
                                                                                  f = (s / u(i,i)) / g

                                                                                  do 450 k = i, m
                                                                                     u(k,j) = u(k,j) + f * u(k,i)
450                                                                                  continue

460                                                                                  do 470 j = i, m
470                                                                                     u(j,i) = u(j,i) / g

                                                                                        goto 490

475                                                                                     do 480 j = i, m
480                                                                                        u(j,i) = 0.0

490                                                                                        u(i,i) = u(i,i) + 1.0
500                                                                                        continue
!!$c       .....diagonalization of the bidiagonal form ....
!!$c       .... for k = n step -1 until 1 do -- .....
510                                                                                        do 700 kk = 1, n
                                                                                              k1 = n - kk
                                                                                              k = k1 + 1
                                                                                              its = 0
!!$       .... test for splitting ....
!!$c               for l = k step -1 until 1 do -- .....
520                                                                                           do 530 ll = 1, k
                                                                                                 l1 = k - ll
                                                                                                 l = l1 + 1
                                                                                                 if(abs(rv1(l)) + anorm .eq. anorm) goto 565
!!$c       .... rv1(i) is always zero. so there is no exit through the
!!$c               bottom of the loop.....
                                                                                                 if(abs(w(l1)) + anorm .eq. anorm) goto 540
530                                                                                              continue
!!$c       .....cancellation of rv1(l) if l greater than 1 .....
540                                                                                              c = 0.0
                                                                                                 s = 1.0

                                                                                                 do 560 i = l, k
                                                                                                    f = s * rv1(i)
                                                                                                    rv1(i) = c* rv1(i)
                                                                                                    if( abs(f) + anorm .eq. anorm) goto 565
                                                                                                    g = w(i)
                                                                                                    h = sqrt(f*f + g*g)
                                                                                                    w(i) = h
                                                                                                    c = g / h
                                                                                                    s = -f / h
                                                                                                    if(.not.matu) goto 560

                                                                                                    do 550 j = 1, m
                                                                                                       y = u(j,l1)
                                                                                                       z = u(j,i)
                                                                                                       u(j,l1) = y * c + z * s
                                                                                                       u(j,i) = - y * s + z * c
550                                                                                                    continue

560                                                                                                    continue
!!$c       .... test for convergence ....
565                                                                                                    z = w(k)
                                                                                                       if(l .eq. k) goto 650
!!$c       .... shift from bottom by 2 by 2 minor ....
                                                                                                       if(its .eq. 30 ) goto 1000
                                                                                                       its = its + 1
                                                                                                       x = w(l)
                                                                                                       y = w(k1)
                                                                                                       g = rv1(k1)
                                                                                                       h = rv1(k)
                                                                                                       f = ((y - z) * ( y + z ) + (g - h) * (g + h))/( 2.0 * h * y)
                                                                                                       g = sqrt ( f*f + 1.0 )
                                                                                                       f = (( x - z) * (x + z) + h * (y / (f + sign(g,f)) - h)) / x
!!$       .... next QR transformation....
                                                                                                       c = 1.0
                                                                                                       s = 1.0

                                                                                                       do 600 i1 = l, k1
                                                                                                          i = i1 + 1
                                                                                                          g = rv1(i)
                                                                                                          y = w(i)
                                                                                                          h = s * g
                                                                                                          g = c * g
                                                                                                          z = sqrt(f*f + h*h)
                                                                                                          rv1(i1) = z
                                                                                                          c = f / z
                                                                                                          s = h / z
                                                                                                          f = x * c + g * s
                                                                                                          g = -x * s + g * c
                                                                                                          h = y * s
                                                                                                          y = y * c
                                                                                                          if(.not.matv) goto 575

                                                                                                          do 570 j = 1, n
                                                                                                             x = v(j,i1)
                                                                                                             z = v(j,i)
                                                                                                             v(j,i1) = x * c + z * s
                                                                                                             v(j,i) = - x * s + z * c
570                                                                                                          continue

575                                                                                                          z = sqrt(f*f + h*h)
                                                                                                             w(i1) = z
!!$c       ....rotation can be arbitrary if z is zero....
                                                                                                             if( z .eq. 0.0) goto 580
                                                                                                             c = f / z
                                                                                                             s = h / z
580                                                                                                          f = c * g + s * y
                                                                                                             x = -s * g + c * y
                                                                                                             if(.not.matu) goto 600

                                                                                                             do 590 j = 1, m
                                                                                                                y = u(j,i1)
                                                                                                                z = u(j,i)
                                                                                                                u(j,i1) = y * c + z * s
                                                                                                                u(j,i) = -y * s + z * c
590                                                                                                             continue

600                                                                                                             continue

                                                                                                                rv1(l) = 0.0
                                                                                                                rv1(k) = f
                                                                                                                w(k) = x
                                                                                                                goto 520
!!$c       .... convergence ...
650                                                                                                             if(z. ge. 0.0) goto 700
!!$c       .... w(k) is made non-negative ....
                                                                                                                w(k) = -z
                                                                                                                if(.not.matv) goto 700

                                                                                                                do 690 j = 1, n
690                                                                                                                v(j,k) = -v(j,k)

700                                                                                                                continue
                                                                                                                   goto 1001
!!$c       .....set error -- no convergence to A
!!$c               singular value after 30 iterations ...
1000                                                                                                               ierr = k
1001                                                                                                               return
                                                                                                                end do
#endif

                                                                                                                
 end Module useful_functions


