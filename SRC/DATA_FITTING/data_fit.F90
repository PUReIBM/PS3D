Program data_fit 

  USE precision  
  USE constants  
  USE fitting_functions
  implicit none   
  INTEGER, PARAMETER :: nimax = 20, njmax = 20, ndatamax = 20, mamax = 10
  INTEGER, parameter :: minunitno = 200, maxunitno = 400 
  INTEGER, parameter ::mpmax = njmax, npmax = mamax
  INTEGER, parameter :: maxfiles = 10

  CHARACTER*80 :: header_ni(nimax) 
  CHARACTER*10 :: xvar, yvar
  REAL(prcn) :: xval(nimax), yval(nimax, njmax), data(nimax, njmax,ndatamax), data2(nimax,njmax), sdev(nimax, njmax), datafit(njmax, mamax), Rsq(nimax)
  INTEGER :: nfiles, ifile , iunit, ni, nj, loopcount, data_nicase(nimax)  
  CHARACTER*200 :: filelistname, ifilename, allfiles(maxfiles)  

  LOGICAL :: nuss_fit
!!$  INTERFACE
!!$     FUNCTION funpoly(x,n)
!!$       USE nrtype
!!$       IMPLICIT NONE
!!$       REAL(DP), INTENT(IN) :: x
!!$       INTEGER(I4B), INTENT(IN) :: n
!!$       REAL(DP), DIMENSION(n) :: funpoly
!!$     END FUNCTION funpoly
!!$  END INTERFACE

  CALL calculate_constants
  CALL fit_data
Contains 
  Subroutine fit_data
    USE general_funcs
    !USE nr , ONly :  svdfit 
    implicit none 
    CHARACTER*200 :: JUNK, tmp_yval, ffilename


    INTEGER :: munit, i, j, n, nj2
    LOGICAL :: EOF, EOZ 
    REAL(prcn) :: chisq, yfit, sse, sst, ymean , tmp, carmandrag, ssdrag, re, phi 
    REAL(prcn), dimension(:,:), allocatable  :: v
    REAL(prcn), dimension(:), allocatable :: w, coeff_a, x, y, sig, xdata, ydata
    INTEGER :: nbasis
    filelistname = "listoffiles.flist" 
    
    !nuss_fit = .TRUE.
    nuss_fit = .FALSE.

    munit=getnewunit(minunitno, maxunitno) 
    CALL openfile(munit, filename=filelistname, formfile="formatted", statusfile = "old")

    READ(munit, '(i4)') nfiles
    IF(nfiles.gt.maxfiles) then 
       WRITE(*,'(A,/,A,/,A)') 'number of files greater than maxfiles used for arrays allocation:', 'Re Complile with higher maxfiles', 'STOPPING THE DATA FIT FOR NOW'
       STOP
    ELSE

       WRITE(*,'(A,i4)') 'NUMBER OF FILES THAT WILL BE OPENED  = ', nfiles
    end IF

    FILES: DO ifile = 1, nfiles
       READ(munit,'(A)') ifilename
       ifilename = TRIM(ifilename) 
       WRITE(*,'(A,A)') 'FILE GOING TO BE OPENED IS:  ', ifilename
       allfiles(ifile) = ifilename


       iunit = getnewunit(minunitno, maxunitno)  

       CALL openfile(iunit, filename = TRIM(ifilename), formfile = 'formatted', statusfile = 'old')
       xvar = ""
       yvar = ""
       read(iunit,*) JUNK

       READ(iunit,'(10x,A)') xvar
       READ(iunit,'(10x,A)') yvar

       write(*,'(3(2x,A))') 'xvar and yvar names: ', TRIM(xvar), TRIM(yvar)
       
       EOF = .false.
       EOZ = .false.
       ni = 0
       
       DO while(.NOT.EOF) 
          
          nj = 0

          EOF = check_eof(iunit)
          IF(EOF) exit 

          EOZ = .false. 
          ni = ni + 1

          DATA(ni,:,:) = zero 
          loopcount = 0
          DO while(.not.EOZ)
             EOF = check_eof(iunit) 
             EOZ = check_eoz(iunit)

             IF(EOF.or.EOZ) exit 
             loopcount = loopcount + 1
             nj = nj + 1

             if(loopcount.eq.1) then 
                read(iunit,'(A)') header_ni(ni) 
                !READ(iunit,'(8x,A)') tmp_yval 
                READ(iunit,'(8x,g17.8)') xval(ni) 
                !xval(ni)

                write(*,'(A,1x,A,1x, g17.8)') 'zone and xvar name = ', header_ni(ni) , xval(ni)
             end if

             read(iunit,*) yval(ni,nj), data(ni,nj,1), data(ni,nj,2), data(ni,nj,3), data(ni,nj,4), data(ni,nj,5)
             
             if(.not.NUSS_FIT)data(ni,nj,4) = data(ni,nj,5) 

             !WRITE(*,'(6(2x,g17.8))')yval(ni,nj), data(ni,nj, 1:5) 
             !EOF = .true.
          end DO
          data_nicase(ni) = nj 
          WRITE(*,'(A,2X,i4,2X,A,2X,i4)') 'NUMBER OF DATA PTS in ZONE', ni, 'EQUAL TO: ', data_nicase(ni)


       end DO


       WRITE(*,'(A,i4)') 'NUMBER OF ZONES IN THIS FILE  = ', ni !data_nicase(nj)

       CLOSE(iunit,status="keep")
       CALL DATA_FIT_2D
       goto 10000

       !Now perform the data fitting 
       ffilename = TRIM(ifilename)//'_fitdata'
       CALL openfile(iunit, filename = TRIM(ffilename), formfile = 'formatted', statusfile = 'unknown')


       do i = 1, ni  
          nj = data_nicase(i) 
          sdev(i, 1:nj) = one
!!$          tmp = data(i,4,3)
!!$          data(i,4,3) = data(i,nj,3)
!!$          data(i,nj,3) = tmp
!!$
!!$          tmp = yval(i,4)
!!$          yval(i,4) = yval(i,nj)
!!$          yval(i,nj) = tmp
          nbasis = 4
          WRITE(*,*) 'xval = ', xval(i) , nj
          nj = nj !- 1 
          nj2 = nj! +1
          ALLOCATE(v(nbasis,nbasis), W(nbasis), coeff_a(nbasis), x(nj), y(nj), sig(nj), xdata(nj2), ydata(nj2))
          y = zero 
          ydata = zero 
          x(1:nj) = (yval(i,1:nj))


          coeff_a = zero
          do j = 1, nj 
             !carmandrag = 10.d0*x(j)/((one-x(j))**3.d0)
             !carmandrag = ((one-x(j))**3.d0)
             !y(j) = (data(i,j,1)  -one)*carmandrag
             !WRITE(*,*) data(i,j,1), carmandrag,data(i,j,1)  -one-  carmandrag

             phi = xval(i) 
             re = x(j) 
             ssdrag = zero
             y(j) = (data(i,j,4))

          end do

          !x(1:nj) = log(x(1:nj))
          xdata(1:nj2) = x(1:nj) !yval(i,1:nj2)
          !do j = 1, nj 
          !y(j) = 3.2*x(j)/((one-x(j))**3.d0)+ 1.7*(x(j)-one) + 1.2*SIN((x(j))**(0.25d0))
          !end do
          !y(1) = zero 
          !ydata(1) = zero

          sig(1:nj) = one 
          v = zero
          w = zero 
          ymean = SUM(y(1:nj))/real(nj,prcn)
          write(*,'(A,i4,/,"yvals = ", 10(2x,g17.8))')'nj = ', nj, x(1:nj)
          write(*,'("data  = ", 10(2x,g17.8))') y(1:nj)
          chisq = zero 

          !CALL svdfit(x,y,sig, coeff_a, v, w,chisq,funstokes)


          WRITE(*,*) 'DATA FITTING ZONE', i
          WRITE(*,'(A,10(1x,g17.8))') 'coeffs = ', coeff_a(1:nbasis)

          write(iunit,'(A)') header_ni(i) 
          do j=1,nj2
             !datafit(j,:)=funstokes (xdata(j),nbasis)
             yfit = zero 
             sse = zero 
             sst = zero 
             do n = 1, nbasis
                yfit = yfit + datafit(j,n)*coeff_a(n)
             end do
             sse = sse + (yfit-ydata(j))**two
             !SSE: sum of the squares residual 
             sst = sst + (ydata(j)-ymean)**two
             !SST: TOTAL SUM OF THE SQUARES
             WRITE(iunit,'(3(2x,g17.8))') xdata(j), ydata(j), yfit
          end do
          Rsq(i) = one - SSE/SST
          WRITE(*,'("goodness of the fit: CHISQ AND RSQ= ", 2(2x,g17.8))') chisq, Rsq(i) 
          WRITE(iunit,'("#coeffs =",10(1x,g17.8))') coeff_a(1:nbasis)

          WRITE(iunit,'("#goodness of the fit: CHISQ AND RSQ= ", 2(2x,g17.8))') chisq, Rsq(i) 

          DEALLOCATE(v, W, coeff_a, x, y, sig,xdata, ydata)
       end do
10000  continue
       close(iunit,status = "keep")
    end DO FILES

    !open(unit = munit, file=filename, form="formatted") 

    !open(unit = funit, file=filename, form="formatted") 

    !READ(funit,*) "junk" 
    !READ(funit) 

  end Subroutine fit_data



  SUBROUTINE DATA_FIT_2D
    USE general_funcs
    !USE nr , ONly :  svdfit 

    implicit none

    INTEGER, PARAMETER :: nmax = 1000, dimn_ind = 2, dimn_dep = 1
    REAL(prcn) :: indvar(nmax,dimn_ind), depvar(nmax, dimn_dep), sig(nmax), yfitma(nmax,mamax), yfit(nmax) , ymean , ydata(nmax,dimn_dep), xdata(nmax,dimn_dep), ydataij(nmax,nmax, dimn_dep), LHS(nmax, mamax), RHS(nmax), SOL(mamax)
    INTEGER :: i,j,k, count_data, nbasis, n, iunit2, ii, jj 
    CHARACTER*200 :: ffilename, formfile, statusfile 
    REAL(prcn) :: ssdrag, re, phi, chisq, rsq, sse, sst, reub, denvhf, tmp, betanot, ress, den , ssnuss, prub 
    REAL(prcn), dimension(:), allocatable :: w, coeff_a
    REAL(prcn), dimension(:, :), allocatable ::v
    LOGICAL :: refile, volfile 
    nbasis = 6
    ALLOCATE(v(nbasis,nbasis), W(nbasis), coeff_a(nbasis))

    ffilename = TRIM(ifilename)//'_PREFIT.dat'
    CALL openfile(iunit, filename = TRIM(ffilename), formfile = 'formatted', statusfile = 'unknown')

    ffilename = TRIM(ifilename)//'_ULTA_PREFIT.dat'
    refile = .false.
    volfile = .false.
    ymean = zero 
    count_data = 0
    do i = 1, ni  
       nj = data_nicase(i) 
       write(iunit,'(A)') header_ni(i) 
       do j = 1, nj 
          count_data = count_data + 1
          if(count_data.gt.nmax) then 
             WRITE(*,'(A/A)') 'NO. OF DATA PTS GT THAN ARRAY DIMENSIONS IN DATA_FIT_2D', ' STOPPING'
             STOP
          end if

          sdev(i, 1:nj) = one

          if(ifile.eq.1) then
             !for files that are versus RE
             indvar(count_data, 1)  = xval(i)          

             indvar(count_data, 2) = yval(i,j)

             re = indvar(count_data,2) 
             phi = indvar(count_data,1) 
             refile =  .true.
             xdata(j,1) = re
          ELSE
             !for files that are vs VOL_FRAC1
             indvar(count_data, 1)  = yval(i,j)          

             indvar(count_data, 2) = xval(i)
             volfile = .true.
             re = indvar(count_data,2) 
             phi = indvar(count_data,1) 
             xdata(j,1) = phi
          end if

          !carmandrag = 10.d0*x(j)/((one-x(j))**3.d0)
          !carmandrag = ((one-x(j))**3.d0)
          !y(j) = (data(i,j,1)  -one)*carmandrag
          !WRITE(*,*) data(i,j,1), carmandrag,data(i,j,1)  -one-  carmandrag
          ssdrag = zero
          !WRITE(*,*) 'VOLFRAC, RE = ', phi, re
          reub = re!/(one-phi)**1.d0
          prub = 0.7d0 
          den = (one)/((one-phi))
          ress = re*den
          ssdrag = (single_sphere_iner_drag(reub))
          ssnuss = single_sphere_iner_nuss(reub,prub) 
          ssdrag = ssdrag!/((one-phi)**(3.d0))
          betanot = 0.313d0*log(ress)/log(one+ress)
          !betanot = betanot + one/((one-phi)**3.d0) - one + phi**(1.d0/3.d0)/((one-phi)**3.d0)  
          !betanot = betanot/((one-phi)**3.d0)
          betanot = betanot!*(one)/((one-phi)**3.d0)
          !betanot = betanot*(10.d0**(2.2d0*(one-phi)/(re)))!*((one-phi)**4.d0)
          !betanot = betanot*re!*exp(7.d0*phi-1.d0)
          tmp = one + 0.06*re - 8.642E-05*re*re
          !WRITE(*,*) 'tmp = '
          if(phi.eq.zero) then 
             !WRITE(*,*) 'HERE'
             !READ(*,*)
             depvar(count_data,1) = one*tmp !ssdrag
             
             depvar(count_data,1) = -log(depvar(count_data,1))/log(one+re) 
             depvar(count_data,1) = depvar(count_data,1)/betanot
             depvar(count_data,1) = depvar(count_data,1) !- one 
             !depvar(count_data,1) = log(data(i,j,1))!/(one+ssdrag)!/(re)!/((one-phi)**2.d0)!/(one+ssdrag)

             !For Stokes case
             
             !depvar(count_data,1) =  log(single_sphere_iner_drag(RE)+ one) 

             if(NUSS_FIT)then
                depvar(count_data,1) =  zero!(single_sphere_iner_nuss(RE,0.7d0)+two)
             else
                depvar(count_data,1) =  zero!(single_sphere_iner_drag(RE)+one) !Trial by Sudheer
             end if
             
             
             !depvar(count_data,1) =  log(single_sphere_iner_nuss(RE,0.7d0)+two) - log(two)
             !depvar(count_data,1) =  (single_sphere_iner_nuss(RE,0.7d0)+two)!/NUSTOK(phi)
             !depvar(count_data,1) = -2.d0*phi -0.03d0*re*re  -1.d0*re -2.d0*phi*re
             ydata(j,1) = depvar(count_data,1)
             !          else if(Re.eq.0.01d0) then
             !             depvar(count_data,1) =  log(NUSTOK(phi))-log(two)
             !             WRITE(*,*)' RE = ', re
             !READ(*,*)
          ELSE
             
             !denvhf = one + (10.d0**(2.d0*phi))*tmp
             !denvhf = 24.d0*((one-phi)**3.d0)*denvhf
             !denvhf = denvhf/(0.313d0*re)
             !depvar(count_data,1) = log(data(i,j,4)) 
             if(NUSS_FIT)then
                depvar(count_data,1) = (data(i,j,4))-(single_sphere_iner_nuss(RE,0.7d0)+two)/(1-phi)!**3.d0
             else
                depvar(count_data,1) = (data(i,j,4))-(single_sphere_iner_drag(RE)+one)/(1-phi)**3.d0 
             end if
             !depvar(count_data,1) = -2.d0*phi -0.03d0*re*re  -1.d0*re -2.d0*phi*re
             !depvar(count_data,1) = (data(i,j,4)) - (two)!NUSTOK(phi)
             
             !depvar(count_data,1) = depvar(count_data,1)/((one-phi)**1.d0)!/single_sphere_iner_nuss(RE,0.7d0) 
             !depvar(count_data,1) = tmp*data(i,j,2)*((one-phi)**3.d0)

             !depvar(count_data,1) = log(depvar(count_data,1)/ssdrag)!*((one-phi)**3.d0)

             !depvar(count_data,1) = -depvar(count_data,1)/log(one+re) 

             !depvar(count_data,1) = depvar(count_data,1)/betanot - one 


             
             !depvar(count_data,1) = log(data(i,j,1)*((one-phi)**0.d0))-log(one)!/(one+ssdrag)!/(re)!/((one-phi)**2.d0)!/(one+ssdrag)
             !depvar(count_data,1) = depvar(count_data,1)*(exp((4.d0*phi)))
             !depvar(count_data,1) = depvar(count_data,1)*Re**0.2d0
             
             !depvar(count_data,1) = data(i,j,4)
             ydata(j,1) = MAX(depvar(count_data,1),zero)
             !indvar(j,1) = log(indvar(j,1))
          end if

          !depvar(count_data,1) = phi*phi - phi*re + re 
          !WRITE(*,*)'DEP  =', phi, re, depvar(count_data,1)
          ymean = ymean + exp(depvar(count_data,1) )
          !WRITE(*,'(A,2(1x,g17.8), A, 4(2x,g17.8))') 'Y(J) = for RE and PHI: ', RE, phi,' = ', ydata(j,1)!,  depvar(count_data,1), ssdrag

          !LHS(count_data,1:nbasis)=funvol (indvar(count_data,1:2),nbasis)

          if(NUSS_FIT)then
             LHS(count_data,1:nbasis)=funnuss (indvar(count_data,1:2),nbasis)
          else
             LHS(count_data,1:nbasis)=funstokes (indvar(count_data,1:2),nbasis)
          end if
          
          !WRITE(*,*) 'phi, re, LHS  =', phi, re, lhs(count_data,1:nbasis)
          RHS(count_data)=MAX(depvar(count_data,1),zero)
          !WRITE(*,'(A,10(2x,g17.8))') 'RE,PHI ',indvar(count_data,1:2),LHS(count_data,1:nbasis), RHS(count_data)
          
          !WRITE(iunit,'(10(2x,g17.8))') xdata(j,1), data(i,j,3),data(i,j,4),data(i,j,3)/(two+ssnuss), data(i,j,4)/ssnuss
          WRITE(iunit,'(10(2x,g17.8))') xdata(j,1), ydata(j,1)!, single_sphere_iner_nuss(RE,0.7d0)  ! - phi**0.5d0!,data(i,j,4),data(i,j,3)/(two+ssnuss), data(i,j,4)/ssnuss
          !WRITE(*,'(10(2x,g17.8))') xdata(j,1), ydata(j,1), ydata(j,1) - ydata(1,1) ! - phi**0.5d0!,data(i,j,4),data(i,j,3)/(two+ssnuss), data(i,j,4)/ssnuss
       end do

       v = zero 
       w = zero 
       coeff_a = zero 
       sig = one 
!!$
!!$       !ydata(1:nj,1) = log(ydata(1:nj,1))
!!$       CALL svdfit(xdata(1:nj,1:1),ydata(1:nj,1),sig(1:nj), coeff_a, v, w,chisq,funre)
!!$       WRITE(*,'(A,10(1x,g17.8))') 'coeffs = ', coeff_a(1:nbasis)

!!$       do j = 1,nj 
!!$          yfitma(j,:)=funre (xdata(j,1:1),nbasis)
!!$          yfit(j) = zero 
!!$          do n = 1, nbasis
!!$             yfit(j) = yfit(j) + yfitma(j,n)*coeff_a(n)
!!$          end do
!!$
!!$
!!$          
!!$          WRITE(iunit,'(10(2x,g17.8))') xdata(j,1) ,ydata(j,1), yfit(j), ABS((yfit(j)-ydata(j,1))/(ydata(j,1)+small_number))*100
!!$       end do
!!$       WRITE(*,'("goodness of the fit: CHISQ= ", 2(2x,g17.8))') chisq
!!$
!!$       WRITE(iunit,'("#goodness of the fit: CHISQ = ", 2(2x,g17.8))') chisq
!!$       
!!$       WRITE(iunit,'("#coeffs",10(2x,g17.8))')  coeff_a(1:nbasis)

!!$
!!$       

    end do
    
    CALL  xwnnls(LHS(1:count_data,1:nbasis),RHS(1:count_data), SOL(1:nbasis),count_data, nbasis)
    
    WRITE(*,'(A,2x,i4)') 'COUNT_DATA = ', COUNT_DATA
    v = zero 
    w = zero 
    coeff_a = zero
    chisq = zero 
    coeff_a(1:nbasis) = SOL(1:nbasis)
    sig = one 
    ymean  = ymean/real(count_data,prcn) 
    !CALL svdfit(indvar(1:count_data,1:2),depvar(1:count_data,1),sig(1:count_data), coeff_a, v, w,chisq,funstokes)
    
    CLOSE(iunit,status="keep")
    ffilename = TRIM(ifilename)//'_POSTFIT.dat'
    CALL openfile(iunit, filename = TRIM(ffilename), formfile = 'formatted', statusfile = 'unknown')
    j = 0 
    do i  = 1, ni 
       nj = data_nicase(i) 
       
       write(iunit,'(A)') header_ni(i) 
       do jj = 1, nj
          j = j + 1
          if(refile) then 
             re = indvar(j,2) 
             phi = indvar(j,1) 
          else
             re = indvar(j,1) 
             phi = indvar(j,2) 
          end if
          
          if(NUSS_FIT)then
             yfitma(j,1:nbasis)=funnuss (indvar(j,1:2),nbasis)
          else
             yfitma(j,1:nbasis)=funstokes (indvar(j,1:2),nbasis)
          end if
          
          yfit(j) = zero 
          sse = zero 
          sst = zero 
          do n = 1, nbasis
             yfit(j) = yfit(j) + yfitma(j,n)*coeff_a(n)
          end do
          if(NUSS_FIT)then
             yfit(j) = yfit(j) + (single_sphere_iner_nuss(RE,0.7d0)+two)/(1-phi)!**3.d0
             depvar(j,1) = depvar(j,1) + (single_sphere_iner_nuss(RE,0.7d0)+two)/(1-phi)!**3.d0
          else
             yfit(j) = yfit(j) + (single_sphere_iner_drag(RE)+one)/(1-phi)**3.d0
             depvar(j,1) = depvar(j,1) + (single_sphere_iner_drag(RE)+one)/(1-phi)**3.d0
          end if
          !if(re.eq.0.01d0) then 
          !   yfit(j)  =NUSTOK(phi)
          !   WRITE(*,*) 'HERE, phi - ', phi
          !   READ(*,*)
          !ELSE 
          !yfit(j) = 2.d0*exp(yfit(j)) 
          !end if
       
          
          !depvar(j,1) = 2.d0*exp(depvar(j,1)) 
          chisq = chisq + (depvar(j,1)-yfit(j))**2.d0
          sse = sse + (yfit(j)-depvar(j,1))**two
          !SSE: sum of the squares residual 
          sst = sst + (depvar(j,1)-ymean)**two
          !SST: TOTAL SUM OF THE SQUARES
          !WRITE(iunit,'(20(2x,g17.8))') indvar(j,1:2), exp(depvar(j,1)), exp(yfit(j)), ABS((exp(depvar(j,1)) -exp( yfit(j)))/(exp(depvar(j,1))+SMALL_NUMBER))*100.d0!, yfitma(j,1:nbasis)
          !WRITE(*,'(20(2x,g17.8))') indvar(j,1:2), depvar(j,1), yfit(j), ABS((depvar(j,1) - yfit(j))/(depvar(j,1)+SMALL_NUMBER))*100.d0!, yfitma(j,1:nbasis)
          WRITE(iunit,'(20(2x,g17.8))') indvar(j,1:2), depvar(j,1), yfit(j), ABS(((depvar(j,1)) -( yfit(j)))/((depvar(j,1))+SMALL_NUMBER))*100.d0!, yfitma(j,1:nbasis)
          !WRITE(iunit,'(20(2x,g17.8))') indvar(j,1:2), MAX(depvar(j,1),zero), yfit(j), ABS((MAX(depvar(j,1),zero) -( yfit(j)))/(MAX(depvar(j,1),zero)+SMALL_NUMBER))*100.d0!, yfitma(j,1:nbasis)
          !WRITE(*,'(20(2x,g17.8))') indvar(j,1:2), depvar(j,1), yfit(j), ABS((depvar(j,1) - yfit(j))/(depvar(j,1)+SMALL_NUMBER))*100.d0!, yfitma(j,1:nbasis)
       end do
    end do

    Rsq = one - SSE/SST
    WRITE(*,'("goodness of the fit: CHISQ AND RSQ= ", 2(2x,g17.8))') chisq, Rsq
    WRITE(iunit,'("#coeffs =",10(1x,g17.8))') coeff_a(1:nbasis)

    WRITE(iunit,'("#goodness of the fit: CHISQ AND RSQ= ", 2(2x,g17.8))') chisq, Rsq
    WRITE(*,'(A,10(1x,g17.8))') 'coeffs = ', coeff_a(1:nbasis)
    CLOSE(iunit, status="keep")

    DEALLOCATE(v, W, coeff_a) 
  end SUBROUTINE DATA_FIT_2D
  SUBROUTINE  xwnnls(vprime,aprime,gamma,mv, nv)
    IMPLICIT NONE
!!$    This subroutine is used to estimate gamma by casting A' = -
    !!\gamma v' into a  linear least squares problem. The approximate
    !! set of equations to be solved are                             
    !!      V'\Gamma = -A' 
!!$  V' => mv X nv,  \Gamma => nv X 1,   A' => mv X 1
!!$ The equations are solved with the non negativity constraints in
    !! the last few rows of the \Gamma vector.

    INTEGER , INTENT(in) :: mv,nv
    INTEGER :: me,lw,liw,L,mdw,mode, row
    REAL(prcn), Intent(out) :: gamma(nv)
    REAL(prcn) :: prgopt(1),upval
    REAL(prcn),INTENT(inout) ::  vprime(mv,nv), aprime(mv)
    REAL(prcn)  :: rnorm
    REAL(prcn) :: vprime_temp(mv,nv), aprime_temp(mv),gamma_temp(nv)
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
    REAL(prcn), DIMENSION(:,:),ALLOCATABLE ::  E
    REAL(prcn), DIMENSION(:),ALLOCATABLE ::  F
    REAL(prcn), DIMENSION(:), ALLOCATABLE :: work
    REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: W


    me = 0 ! For the input to DWNNLS

!!$    Interchange rows and columns corresponding to \gamma_11, 
    !!\gamma_22 and \gamma_33 (since non negativity constraints are
    !! imposed on these values)
    vprime_temp = vprime
!!$    do row = 1, mv
!!$       vprime(row,1) = vprime_temp(row,nv-ndim+1) ! v'_i1 --> v'_i7
!!$       vprime(row,nv-ndim+1) = vprime_temp(row,1) ! v'_i7 --> v'_i1
!!$       vprime(row,1+ndim+1) = vprime_temp(row,nv-ndim+2) ! v'_i5 --> v'_i8
!!$       vprime(row, nv-ndim+2) = vprime_temp(row,1+ndim+1) ! v'_i8 --> v'_i5
!!$    end do
    aprime_temp = aprime
    mdw  = me+mv
    !k = max(ma+mg,nmat)
    lw = me+mv+5*nv
    liw = me+mv+nv
    IF(.NOT.ALLOCATED(iwork)) ALLOCATE(iwork(liw))
    IF(.NOT.ALLOCATED(work)) ALLOCATE(work(mv+5*nv))
    IF(.NOT.ALLOCATED(W)) ALLOCATE(W(mdw,nv+1))
    IF(me.GT.0) THEN 
       W(1:me,1:nv) = E(1:me,1:nv)
       W(1:me,nv+1) = F(1:me)
    ENDIF
    W(me+1:me+mv,1:nv) = vprime(1:mv,1:nv)
    W(me+1:me+mv,nv+1) = aprime(1:mv)
    iwork(1) = lw
    iwork(2) = liw
    prgopt(1) = 1.
    !PRINT*,'MA = ', MV
    !PRINT*,'NA = ', NV
    !L = 3
    !L=nv-ndim
    L = nv
    !PRINT*,'L = ', L
    CALL DWNNLS (W, MDW, ME, MV, NV, L, PRGOPT, gamma_temp, RNORM, MODE, &
         IWORK, WORK)
    gamma(1:nv) = gamma_temp(1:nv)
    !IF(mode.NE.0)   
    !PRINT*,'Mode=',mode

!!$    gamma(1) = gamma_temp(nv-ndim+1)
!!$    gamma(nv-ndim+1) = gamma_temp(1)
!!$    gamma(1+ndim+1) = gamma_temp(nv-ndim+2)
!!$    gamma(nv-ndim+2) = gamma_temp(1+ndim+1)
    do row=1,nv
       PRINT*,'gamma = ' ,gamma(row)
    end do

  end SUBROUTINE xwnnls

    
  LOGICAL FUNCTION check_eoz(funit) 
    Implicit none
    integer, intent(in) :: funit
    CHARACTER*20 :: rstring
    read(funit,'(1x,A)') rstring
    IF(TRIM(rstring).eq."EOZ") then 
       check_eoz = .true. 
    ELSE
       !rewind the file by one record if end of zone not yet reached 
       check_eoz = .false.
       backspace(funit)
    end IF

  end FUNCTION check_eoz

  LOGICAL FUNCTION check_eof(funit) 
    Implicit none
    integer, intent(in) :: funit
    CHARACTER*20 :: rstring
    read(funit,'(1x,A)') rstring
    IF(TRIM(rstring).eq."EOF") then 
       check_eof = .true. 
    ELSE
       !rewind the file by one record if end of file not yet reached 
       check_eof = .false.
       backspace(funit) 
    end IF

  end FUNCTION check_eof

  REAL(prcn) FUNCTION NUSTOK(phi) 
    Implicit none 
    REAL(prcn) :: phi
    REAL(prcn) :: tmp 
    tmp = one/((one-phi)**three)
    
    NUSTOK  = two + tmp*(14.14d0*phi + 3.89d0*phi*phi - 21.23*(phi**1.5d0) + 1.09d0*(phi**(one/three)))
    
    RETURN
  end FUNCTION NUSTOK

  REAL(prcn) FUNCTION FSTOK(phi) 
    Implicit none 
    REAL(prcn) :: phi
    REAL(prcn) :: tmp 
    tmp = (one-phi)
    
    FSTOK  = (10.d0*phi/(tmp**two) + (tmp**two)*(1+1.5*dsqrt(phi)))/tmp
    PRINT*, 'phi =  ', phi, 'F = ', FSTOK
    RETURN
  end FUNCTION FSTOK

  REAL(prcn) FUNCTION single_sphere_iner_drag(rein) 
    Implicit none 
    REAL(prcn) :: rein
    single_sphere_iner_drag = 0.15d0*(REIN**0.687d0)
  end FUNCTION single_sphere_iner_drag

  REAL(prcn) FUNCTION single_sphere_iner_nuss(rein,prin) 
    Implicit none 
    REAL(prcn) :: rein, prin
    REAL(prcn) :: tmp
    tmp = zero  
    if(rein.ge.one) then 
       tmp = rein*prin
       tmp = one + one/tmp
       tmp = (tmp**(one/three))*(rein**0.41d0)
       tmp = tmp*(prin**(one/three))
       tmp = one + tmp 

       
    ELSE
       tmp = 0.6d0*(rein**0.5d0)*(prin**(one/three))
       tmp = tmp + two 
    end if
    
    single_sphere_iner_nuss = tmp - two 
    
  end FUNCTION single_sphere_iner_nuss



end Program data_fit



