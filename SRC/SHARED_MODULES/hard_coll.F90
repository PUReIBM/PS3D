MODULE hard_sphere 
  USE global_data, ONLY : RUN_NAME, gofr_avg, rad_bin, nbins
  USE precision 
  USE constants 
  USE randomno 
  Private 
  INTEGER, PARAMETER :: NDIM=3
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  rx,ry,rz,vx, vy, vz&
       &, coltim, rp, rx_0, ry_0, rz_0, rxt, ryt, rzt, wt, wt1,wt2, ftemp1, ftemp2 
  DOUBLE PRECISION rhod,vol,t, radp1,radp2,b_colimp, rest, ener0
  DOUBLE PRECISION delenergy, ener, t_ener,r_0
  DOUBLE PRECISION ener1, ener2, ener3, mnsq0, kbt,mnfrepth
  DOUBLE PRECISION ener_spec1, ener_spec2,mass_weighted_var,mweighted_var0,enersp1_init,enersp2_init
  DOUBLE PRECISION rmsvel, diff_kth, diff_kth2, cse_gofr
  DOUBLE PRECISION t_numcoll, variance, m_p, lambda, variance1,  variance2
  DOUBLE PRECISION diff_ensk, t_numkth1, t_numkth2, t_numensk
  INTEGER, ALLOCATABLE, DIMENSION(:) :: partnr, ifpon, itag, icolon
  INTEGER    coll,ihist,ncoal,ngraz,igraz !,nbins
  INTEGER    iwriteout,ikivaecoal
  INTEGER      ncoll,  overlap_count, n
  
  Logical bump_in_uplist, overlap_uplist 
  PUBLIC :: perform_hard_sphere
!!$        common / block1 / rx, ry, rz, vx, vy, vz,t, radp1, radp2, rest
!!$        common / block2 / coltim,wt,rp,rhod,vol,frandr,tij,b_colimp
!!$        common / block3 / ener0, delenergy, ener, ncoll
!!$	common / block4 / t_ener, r_0, lambda
!!$	common / block5 / diff_kth, diff_kth2, diff_ensk, cse_gofr
!!$	common / block6 / t_numkth1, t_numkth2, t_numensk
!!$	common / block7 / ener1, ener2, ener3
!!$	common / block8 / t_numcoll, variance, m_p, mnsq0
!!$	common / block9 / rx_0, ry_0, rz_0, rxt, ryt, rzt, kbt,mnfrepth


!!$	common /intblock/ partnr,n,coll,nbins,ihist,ncoal,ngraz,igraz
!!$	common /intblock2/ icolnon,iwriteout,ifpon,ikivaecoal,ibidisperse
!!$        common /intblock3/itag

  ! *******************************************************************
  ! ** THIS FORTRAN CODE IS INTENDED TO ILLUSTRATE POINTS MADE IN    **
  ! ** THE TEXT. TO OUR KNOWLEDGE IT WORKS CORRECTLY. HOWEVER IT IS  **
  ! ** THE RESPONSIBILITY OF THE USER TO TEST IT, IF IT IS USED IN A **
  ! ** RESEARCH APPLICATION.                                         **
  ! *******************************************************************

  ! *******************************************************************
  ! ** FICHE F.10  MOLECULAR DYNAMICS OF HARD SPHERE ATOMS           **
  ! *******************************************************************
CONTAINS 
  SUBROUTINE perform_hard_sphere(Npart,nbody1,nbody2,coeff_e, pvel_var1, pvel_var2, xc, Tratio, radbdy, partdens, MINCOLS, MIS, ibidisperse)
    USE general_funcs
    USE postproc_funcs

    IMPLICIT NONE 
    
    ! *******************************************************************
    ! ** THIS PROGRAM TAKES IN A HARD-SPHERE CONFIGURATION (POSITIONS  **
    ! ** AND VELOCITIES), CHECKS FOR OVERLAPS, AND THEN CONDUCTS A     **
    ! ** MOLECULAR DYNAMICS SIMULATION RUN FOR A SPECIFIED NUMBER OF   **
    ! ** COLLISIONS.  THE PROGRAM IS FAIRLY EFFICIENT, BUT USES NO     **
    ! ** SPECIAL NEIGHBOUR LISTS, SO IS RESTRICTED TO A SMALL NUMBER   **
    ! ** OF PARTICLES (<500).  IT IS ALWAYS ASSUMED THAT COLLISIONS    **
    ! ** CAN BE PREDICTED BY LOOKING AT NEAREST NEIGHBOUR PARTICLES IN **
    ! ** THE MINIMUM IMAGE CONVENTION OF PERIODIC BOUNDARIES.          **
    ! ** THE BOX IS TAKEN TO BE OF UNIT LENGTH.                        **
    ! ** HOWEVER, RESULTS ARE GIVEN IN UNITS WHERE SIGMA=1, KT=1.      **
    ! **                                                               **
    ! ** PRINCIPAL VARIABLES:                                          **
    ! **                                                               **
    ! ** INTEGER N                    NUMBER OF ATOMS                  **
    ! ** REAL    RX(N),RY(N),RZ(N)    ATOM POSITIONS                   **
    ! ** REAL    VX(N),VY(N),VZ(N)    ATOM VELOCITIES                  **
    ! ** REAL    COLTIM(N)            TIME TO NEXT COLLISION           **
    ! ** INTEGER PARTNR(N)            COLLISION PARTNER                **
    ! ** REAL    SIGMA                ATOM DIAMETER                    **
    ! **                                                               **
    ! ** ROUTINES REFERENCED:                                          **
    ! **                                                               **
    ! ** SUBROUTINE READCN ( CNFILE )                                  **
    ! **    READS IN CONFIGURATION                                     **
    ! ** SUBROUTINE CHECK ( SIGMA, OVRLAP, E )                         **
    ! **    CHECKS CONFIGURATION AND CALCULATES ENERGY                 **
    ! ** SUBROUTINE UPLIST ( SIGMA, I )                                **
    ! **    SEEKS COLLISIONS WITH J>I                                  **
    ! ** SUBROUTINE DNLIST ( SIGMA, I )                                **
    ! **    SEEKS COLLISIONS WITH J<I                                  **
    ! ** SUBROUTINE BUMP ( SIGMA, I, J, W )                            **
    ! **    DOES COLLISION DYNAMICS AND CALCULATES COLLISION VIRIAL    **
    ! ** SUBROUTINE WRITCN ( CNFILE )                                  **
    ! **    WRITES OUT CONFIGURATION                                   **
    ! *******************************************************************

    INTEGER , Intent(in) ::     Npart, mincols, MIS, nbody1, nbody2
    !PARAMETER ( N = 108 )
    INTEGER, PARAMETER :: nbin_max = 201, nhbins=20
    LOGICAL, Intent(in) :: ibidisperse
    DOUBLE PRECISION, INTENT(IN) :: coeff_e, pvel_var1, pvel_var2
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(Npart,3) :: xc
    DOUBLE PRECISION, INTENT(INOUT) :: Tratio
    DOUBLE PRECISION, INTENT(IN), DIMENSION(Npart) :: radbdy
    DOUBLE PRECISION, INTENT(IN) :: partdens
    double precision tst, tfinal, vtemp(npart,ndim), umf0(ndim), rest_input
    double precision rsf(ndim,ndim),rsf1(ndim,ndim),rsf2(ndim,ndim), alpha, t_freepath, eneranaly,&
         & conf1, conf2 , &
         & rho_est_avg(nbin_max), gofr(mis, nbin_max), rho_est(mis,&
         & nbin_max), n_mean(mis,nbin_max), n_var(mis,nbin_max), vol_bin(nbin_max)&
         & , nmean_avg(nbin_max), nvar_avg(nbin_max), vtemp1(nbody1,ndim), vtemp2(nbody2,ndim) !gofr_avg(nbin_max), rad_bin(nbin_max)
   
    LOGICAL :: rescaling, cont_simulation, xperiodicc
    integer colln_per_part(npart), MINIMUM_COLS,  MAXIMUM_COLS, nsim,&
         & nrbins, gofunit, phisq_unit, unit1, unit2, idim !nbins, 
    REAL(prcn) ::       TIMBIG, ftemp(npart), hist(nhbins), fmin, fmax,scale_up1,scale_up2
    PARAMETER ( TIMBIG = 1.0E10 )

    !REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
    !REAL        COLTIM(N)
    !INTEGER     PARTNR(N)
    REAL(prcn)        SIGMA
    
    INTEGER     I, J, K
    REAL(prcn)        DENSTY, DIJ, TIJ, RATE
    REAL(prcn)        E, EN, ENKT, W, PVNKT1, ACW, TEMP, TBC
    CHARACTER   TITLE*80, CNFILE*30
    LOGICAL     OVRLAP
    double precision rxi, ryi, rzi, rxij, ryij, rzij
    double precision vxi, vyi, vzi, vxij, vyij, vzij, bij
    
    nrbins = 200 ! for calculation og g(r)
    nbins = 20   ! for computation of <N^2>
    IF(nrbins.GT.nbin_max) THEN
       WRITE(*,*) "INCREASE nbin_max in hard_coll"
       STOP
    end IF
    IF(nbins.GT.nbin_max) THEN
       WRITE(*,*) "INCREASE nbin_max in hard_coll"
       STOP
    end IF
    
!!$    WRITE(*,*) 'N = ', Npart, ' MINCOLS = ', MINCOLS,' Var = ', pvel_var, 'rest =', coeff_e, MAXVAL(radbdy)
    WRITE(*,*) 'N in HARD SPHERE = ', NPART
    N = Npart
    ! *******************************************************************
    ALLOCATE(rx(n),ry(n),rz(n), vx(n),vy(n), vz(n), wt(n),coltim(n)&
         &,rp(n),rx_0(n),ry_0(n),rz_0(n), rxt(n),ryt(n),rzt(n) )
    ALLOCATE( partnr(n))
    if(ibidisperse) ALLOCATE(wt1(nbody1), wt2(nbody2), ftemp1(nbody1), ftemp2(nbody2))
    
    rest_input = coeff_e
    rest  = one !First Run the elastic case 
    
    xc(1:n,1) = xc(1:n,1) - half
    xc(1:n,2) = xc(1:n,2) - half
    xc(1:n,3) = xc(1:n,3) - half
    rescaling = .FALSE.
    xperiodicc = .FALSE.
    gofunit = getnewunit(100,400)
    
    OPEN(gofunit,file= TRIM(RUN_NAME)//'_gof_avg.dat', form = 'formatted')
    
    phisq_unit = getnewunit(100,400)
    OPEN(phisq_unit,file= TRIM(RUN_NAME)//'_phi_avg.dat', form = 'formatted')
    OPEN(unit=99,file= TRIM(RUN_NAME)//'_hard_coll_related.dat')
    CNFILE = TRIM(RUN_NAME)//"_xc_hs.dat"
    lambda = 1.d0

    
    
    r_0 = radbdy(1) 
    alpha = 0.01d0
    rx(1:N) = xc(1:N, 1)
    ry(1:N) = xc(1:N, 2)
    rz(1:N) = xc(1:N, 3)
    
    rp(1:N) = radbdy(1:N)
    
    rhod = partdens
    if(ibidisperse)then
       variance1 = pvel_var1    
       variance2 = pvel_var2    
       PRINT*,'vel variance1 = ', pvel_var1, 'vel variance2 =', pvel_var2
       do j=1,ndim
          umf0(j)=0.0
          do i=1,ndim
             if(i.eq.j)then
                rsf1(i,j)=variance1
                rsf2(i,j)=variance2
             else
                rsf1(i,j)=0.0
                rsf2(i,j)=0.0
             endif
          enddo
       enddo
       CALL jn_dist(vtemp1,nbody1,ndim,umf0,rsf1)
       CALL jn_dist(vtemp2,nbody2,ndim,umf0,rsf2)
       do i=1,n
          if(i.le.nbody1)then
             vx(i) = vtemp1(i,1)
             vy(i) = vtemp1(i,2)
             vz(i) = vtemp1(i,3)
          else
             vx(i) = vtemp2(i-nbody1,1)
             vy(i) = vtemp2(i-nbody1,2)
             vz(i) = vtemp2(i-nbody1,3)
          end if
       enddo
    else

       variance = pvel_var1    
       do j=1,ndim
          umf0(j)=0.0
          do i=1,ndim
             if(i.eq.j)then
                rsf(i,j)=variance
            else
               rsf(i,j)=0.0
            endif
         enddo
      enddo
      CALL jn_dist(vtemp,n,ndim,umf0,rsf)
      
      do i=1,n
         vx(i) = vtemp(i,1)
         vy(i) = vtemp(i,2)
         vz(i) = vtemp(i,3)
      enddo
   end if
    
    if(ibidisperse)then
       do i=1,nbody1
          wt1(i)=1./float(nbody1)
       end do
       do i=1,nbody2
          wt2(i)=1./float(nbody2)
       end do
       
       unit1 = getnewunit(100, 400)
       OPEN(unit1, FILE= TRIM(RUN_NAME)//'_velpdf1_initial.dat', status='replace')
       unit2 = getnewunit(100, 400)
       OPEN(unit2, FILE= TRIM(RUN_NAME)//'_velpdf2_initial.dat', status='replace')
       do idim=1,ndim
          if(idim.eq.1) then 
             ftemp1(1:nbody1) = vx(1:nbody1)
             ftemp2(1:nbody2) = vx(nbody1+1:n)
          elseif(idim.eq.2) then 
             ftemp1(1:nbody1) = vy(1:nbody1)
             ftemp2(1:nbody2) = vy(nbody1+1:n)
          elseif(idim.eq.3) then 
             ftemp1(1:nbody1) = vz(1:nbody1)
             ftemp2(1:nbody2) = vz(nbody1+1:n)
          end if
      
          CALL histogram(ftemp1,wt1,nbody1,nhbins,fmin,fmax,hist(1:nhbins))
          CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,unit1,real(idim,prcn),1.d0)
          CALL histogram(ftemp2,wt2,nbody2,nhbins,fmin,fmax,hist(1:nhbins))
          CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,unit2,real(idim,prcn),1.d0)
       end do
       close(unit1, status='keep')
       close(unit2, status='keep')
    else
       do i=1,n
          wt(i)=1./float(n)
          colln_per_part(i) = 0
       enddo

   end if
   
    vol = 1.0
    sigma = 2.*r_0
    densty = float(n)*sigma**3
    
!!$    write(*,'('' sigma calculated       '',f15.5)')sigma 
!!$    write(*,'('' reduced density is     '',f15.5)') densty

!!$    PRINT*,'nbfore jn_dist = ', n, size(vtemp,1), size(vtemp,2), rsf, umf0 

    !create weights for histogram calculation

    
    ! ** CHECK FOR PARTICLE OVERLAPS **
    ! ** CALCULATE ENERGY            **

    !CALL CHECK ( SIGMA, OVRLAP, E )

    !IF ( OVRLAP ) STOP 'PARTICLE OVERLAP IN INITIAL CONFIGURATION'
    
    call computeener(nbody1, nbody2)          !(ener2,ener3)
    ener0 = ener
    eneranaly = ener
    if(ibidisperse)then
       enersp1_init = ener_spec1
       enersp2_init = ener_spec2
       mweighted_var0 = mass_weighted_var
    end if
    
!!$    write(*,*)'Done calculating initial energy...'
!!$    write(*,*)'Total energy initial time ..',ener 
    
    !write(*,*)'Energy in 2 components...',ener2,ener3

!!$    open(unit1, FILE='gran_temp.dat', status='replace')
!!$    write(unit1,'(5(E20.10,1x))')0.0, mass_weighted_var/mweighted_var0, ener_spec1/enersp1_init, ener_spec2/enersp2_init, ener_spec2/ener_spec1
!!$    write(16,*)'0.0', ener, '0.0', '0.0'
    
    
    densty = dble(n)*(2.d0*rp(1))**3 
    
    cse_gofr =  (1.-alpha/2.)/(1.-alpha)**3. 
    
    temp = 2./3.*ener0 



    write(*,*)'Variance given to Gaussian random generator.... ',variance  
    write(*,*)'Energy from variance ........................' ,3./2.*variance 
    write(*,*)'Compare with energy from calculation.......',ener0 

!!!!!!!mass of sphere 

    m_p = 1.                  !4./3.*pi*rp(1)**3*rhod 
    write(*,*)'Radius of particle.....',r_0, rp(1) 

    kbt = temp                !!!! m_p/3.d0*variance 

    rmsvel = dsqrt(variance) 

    t_ener = 4.*dble(n)*(2.*r_0)**2.*cse_gofr*sqrt(pi*kbt/m_p) 

    write(99,*)'--------------------------------------------------' 
    write(99,*)'# particles..',n,' radius...',r_0,' Pi..',pi 

    !c     calculate various initial parameters 

    write(99,*)'Mass of sphere...........................',m_p 
    write(99,*)'--------------------------------------------------' 
    write(99,*)'Volume fraction...',alpha,'Density..',densty 
    write(99,*)'--------------------------------------------------' 

    write(99,'((a),2(F20.10,1x))')'Enskog collision frequency...',t_ener 
    write(99,'((a),1(F20.10,1x))')'k_B T ...................... ',kbt 
    write(99,'((a),1(F20.10,1x))')'Velocity variance............',3.* kbt/m_p 

    write(99,'((a),2(F20.10,1x))')'RMS velocity............... ..',sqrt(3.*kbt/m_p) 

    diff_kth  = 1./3.*rmsvel*1./(pi*(2.*r_0)**2.*dble(n)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!$    write(99,'((a),2(F20.10,1x))')'Diffusion coefficient.. .(K Theory)',diff_kth 

    diff_kth2 = 3./8.*1./(dble(n)*(2.*r_0)**2)*sqrt(kbt/(m_p*pi)) !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!$    write(99,'((a),2(F20.10,1x))')'Diffusion coefficient...(K Theory - 2)', diff_kth2 

    !nskog diffusion coefficient for dense regimes (alpha < 0.4) 
    !reference: Mol. Phys. 2003. Vol. 101. No.3. 469-482 
    !Sigurgeirsson and Heyes 

!!$    write(99,'((a),2(F20.10,1x))')'Carnahan Starling g(r).................', cse_gofr 

    diff_ensk = 1.01896/cse_gofr*diff_kth2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!$    write(99,'((a),2(F20.10,1x))')'Diffusion coefficient...(Enskog Theory)', diff_ensk 

    !calculate timescale for scaling number density profile 

    mnfrepth = sqrt(kbt/m_p)*2./t_ener
    t_freepath = mnfrepth**2/diff_kth2

    lambda = lambda/(2.*pi) 

    t_numkth1 = lambda**2/diff_kth 
    t_numkth2 = lambda**2/diff_kth2 
    t_numensk = lambda**2/diff_ensk 
    t_numcoll = 2.d0*lambda**2*t_ener/(variance*cse_gofr) 
!!$
!!$    write(99,'((a),1(F20.10,1x))') 'Numdensity timescale.. kth1.........',t_numkth1 
!!$    write(99,'((a),1(F20.10,1x))')'Numdensity timescale...kth2.........',t_numkth2 
!!$    write(99,'((a),1(F20.10,1x))') 'Numdensity timescale...ensk.........',t_numensk 
!!$    write(99,'((a),1(F20.10,1x))') 'Numdensity timescale...tnumcol......',t_numcoll 
!!$
    
    close(99)
    
    temp = 2.0 * ener / 3.0
    enkt = ener / temp
    !write(*,'('' initial e/nkt        '',f15.5)') enkt

    !** set up initial collision lists coltim and partnr **

    do i = 1, n
       coltim(i) = timbig
       partnr(i) = n
    end do

    open(unit=1000,file=cnfile,form='formatted',status='unknown')
    
    write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
         &    ' "UX" '
    WRITE(1000,*) 'ZONE'
    DO i=1,n
       WRITE(1000,'(4(2x,f12.8))') rx(i), ry(i), rz(i), rp(i)
    enddo
    
    
    do  i = 1, n
       call uplist ( sigma, i )
    end do


    
!call writcn (cnfile)

    !     *******************************************************************
    !     ** main loop begins                                              **
    !     *******************************************************************
    !     ** zero virial accumulator **

    t = 0.0

    !call calgofr       !calculate initial g(r)	

    i = 1
    j = 2 

    ! *******************************************************************
    ! ** MAIN LOOP BEGINS                                              **
    ! *******************************************************************
    MINIMUM_COLS = 0
    MAXIMUM_COLS = 0 
    do while (MINIMUM_COLS.LT.MINCOLS)
       coll = coll + 1

       tij = timbig

       do k = 1, n
          if ( coltim(k) .lt. tij ) then
             tij = coltim(k)
             i   = k
          endif
       enddo

       j = partnr(i)          !i is surely not turned off, nor is j
       
       t = t + tij

       colln_per_part(i) = colln_per_part(i) + 1
       colln_per_part(j) = colln_per_part(j) + 1

       !    ** MOVE PARTICLES FORWARD BY TIME TIJ **
       !    ** AND REDUCE COLLISION TIMES         **
       !    ** APPLY PERIODIC BOUNDARIES          **


       do k = 1, n
          coltim(k) = coltim(k) - tij

          rx(k) = rx(k) + vx(k) * tij
          ry(k) = ry(k) + vy(k) * tij
          rz(k) = rz(k) + vz(k) * tij

          rx(k) = rx(k) - anint ( rx(k) )
          ry(k) = ry(k) - anint ( ry(k) )
          rz(k) = rz(k) - anint ( rz(k) )
       ENDDO

!!$       call computeener(nbody1, nbody2)
!!$       call computeener   !(ener2,ener3)
       
       !     ** compute collision dynamics **
       call bump(sigma,i,j,w) 
       
       
       do k = 1, n
          call uplist ( sigma, k )
       enddo
       
       if(mod(coll,1000).eq.0)write(*,*)'NCOLL DONE..',coll&
            &,colln_per_part(i),colln_per_part(j),t,ener 
       
!!$       if(mod(coll,100).eq.0)then
       call computeener(nbody1, nbody2)
!!$          call computeener   !(ener2,ener3)
!!$       write(unit1,'(5(E20.10,1x))')t, mass_weighted_var/mweighted_var0, ener_spec1/enersp1_init, ener_spec2/enersp2_init, ener_spec1/ener_spec2       
       
       !write(16,'5(E20.10,1x)')t, ener, dble(coll), acw
       !write(49,'7(E20.10,1x)')t,t*t_ener,ener/ener0,float(coll)&
       !     &/t,ener1,ener2,ener3 
    
       MINIMUM_COLS = MINVAL( colln_per_part(1:N))
       MAXIMUM_COLS = MAXVAL( colln_per_part(1:N))
    ENDDO
    !     *******************************************************************
    !     ** main loop ends.                                               **
    !     *******************************************************************

    rx_0(1:n)  = rx(1:n)
    ry_0(1:n)  = ry(1:n)
    rz_0(1:n)  = rz(1:n)

    rest = rest_input !Now run with the real coeff of resitution


    WRITE(1000,*) 'ZONE'
    !write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
    !     &    ' "UX" '
    DO i=1,n
       WRITE(1000,'(4(2x,f12.8))') rx(i), ry(i), rz(i), rp(i)
    enddo
    
   
    nsim  = 0 
    WRITE(*,*) 'DONE WITH ELASTIC CASE: NOW GOING TO REAL CASE'

    overlap_count = 0 
2001 continue 
    
    
    WRITE(*,*) 'NSIM = ', NSIM + 1
    WRITE(*,*) 'coeff. of rest. = ', rest 
    


    do i = 1, n
       coltim(i) = timbig
       partnr(i) = n
       colln_per_part(i)  = 0 !reset the number of collisions for
       ! each realization 
    end do

    !Reset the particle positions to the one obtained from elastic case
    rx(1:n)  = rx_0(1:n)
    ry(1:n)  = ry_0(1:n)
    rz(1:n)  = rz_0(1:n)

!!$    !Renew the velocity distribution 
!!$    CALL jn_dist(vtemp,n,ndim,umf0,rsf)
!!$    
!!$    do i=1,n
!!$       vx(i) = vtemp(i,1)
!!$       vy(i) = vtemp(i,2)
!!$       vz(i) = vtemp(i,3)
!!$    enddo
    
    
    do  i = 1, n
       call uplist ( sigma, i )
    end do
    
!!$    t = 0.0
    
    !call calgofr       !calculate initial g(r)	
    
    i = 1
    j = 2 
    
    MINIMUM_COLS = 0
    MAXIMUM_COLS = 0 
    coll = 0 
    cont_simulation = .TRUE.
!!$    do while (MINIMUM_COLS.LT.MINCOLS)
!!$    do while(mass_weighted_var/mweighted_var0 .gt. 1E-10)
    do while(cont_simulation)
       coll = coll + 1

       tij = timbig

       do k = 1, n
          if ( coltim(k) .lt. tij ) then
             tij = coltim(k)
             i   = k
          endif
       enddo

       j = partnr(i)          !i is surely not turned off, nor is j
       
       t = t + tij

       colln_per_part(i) = colln_per_part(i) + 1
       colln_per_part(j) = colln_per_part(j) + 1

       !    ** MOVE PARTICLES FORWARD BY TIME TIJ **
       !    ** AND REDUCE COLLISION TIMES         **
       !    ** APPLY PERIODIC BOUNDARIES          **


       do k = 1, n
          coltim(k) = coltim(k) - tij

          rx(k) = rx(k) + vx(k) * tij
          ry(k) = ry(k) + vy(k) * tij
          rz(k) = rz(k) + vz(k) * tij

          rx(k) = rx(k) - anint ( rx(k) )
          ry(k) = ry(k) - anint ( ry(k) )
          rz(k) = rz(k) - anint ( rz(k) )
       ENDDO

!!$       call computeener(nbody1, nbody2)
!!$       call computeener   !(ener2,ener3)
       
       !     ** compute collision dynamics **
       call bump(sigma,i,j,w) 
       
       
       do k = 1, n
          overlap_uplist = .false.
          call uplist ( sigma, k )
          IF(overlap_uplist) then 
             overlap_count = overlap_count + 1
             WRITE(*,*) 'OVERLAP TRUE.. REDOING THIS REALIZATION'
             WRITE(*,*) 'No. of REDOS SO FAR = ', overlap_count 
             goto 2001
          end IF
          
       enddo

       if(mod(coll,1000).eq.0)write(*,*)'NCOLL DONE..',coll&
            &,colln_per_part(i),colln_per_part(j),t,ener 

!!$       if(mod(coll,100).eq.0)then
       call computeener(nbody1, nbody2)
!!$          call computeener   !(ener2,ener3)
!!$       write(unit1,'(5(E20.10,1x))')t, mass_weighted_var/mweighted_var0, ener_spec1/enersp1_init, ener_spec2/enersp2_init, ener_spec1/ener_spec2       
       
       !write(16,'5(E20.10,1x)')t, ener, dble(coll), acw
       !write(49,'7(E20.10,1x)')t,t*t_ener,ener/ener0,float(coll)&
       !     &/t,ener1,ener2,ener3 
       
       MINIMUM_COLS = MINVAL( colln_per_part(1:N))
       MAXIMUM_COLS = MAXVAL( colln_per_part(1:N))
       if(rest.eq.one) then
          cont_simulation = (MINIMUM_COLS.LT.MINCOLS)

       else
          cont_simulation = (mass_weighted_var/mweighted_var0 .gt. 1E-10)
       end if
       
    ENDDO !do while 

    xc(1:N, 1) = rx(1:N) + half 
    xc(1:N, 2) = ry(1:N) + half
    xc(1:N, 3) = rz(1:N) + half 
!!$
!!$    velbdy(1:N, 1) = vx(1:N) 
!!$    velbdy(1:N, 2) = vy(1:N) 
!!$    velbdy(1:N, 3) = vz(1:N) 

!!$    if(rest.ne.one.and.ibidisperse)then

!!$    scale_up1 = dsqrt(enersp1_init/ener_spec1)*0.5
!!$    scale_up2 = dsqrt(enersp2_init/ener_spec2)*0.5
!!$       scale_up1 = dsqrt(mweighted_var0/mass_weighted_var)
!!$       scale_up2 = dsqrt(mweighted_var0/mass_weighted_var)
!!$    PRINT*, 'SCALING UP THE VELOCITIES BY', SCALE_UP1, SCALE_UP2
!!$
!!$    velbdy(1:nbody1,1)= velbdy(1:nbody1,1)*scale_up1
!!$    velbdy(nbody1+1:N,1)= velbdy(nbody1+1:N,1)*scale_up2
!!$    velbdy(1:nbody1,2)= velbdy(1:nbody1,2)*scale_up1
!!$    velbdy(nbody1+1:N,2)= velbdy(nbody1+1:N,2)*scale_up2
!!$    velbdy(1:nbody1,3)= velbdy(1:nbody1,3)*scale_up1
!!$    velbdy(nbody1+1:N,3)= velbdy(nbody1+1:N,3)*scale_up2
!!$  endif 
    if(ibidisperse)then
       Tratio = ener_spec2/ener_spec1
       PRINT*,'TRATIO=',Tratio 
    else
       Tratio = zero
    end if
    
    nsim  = nsim + 1
    
    CALL calc_gofr(n, xc(1:n,1:3), radbdy(1:n), 1,1,xperiodicc, nrbins, ndim, rescaling, gofr(nsim,1:nrbins), rho_est(nsim,1:nrbins), rad_bin(1:nrbins)  )
    
    CALL calc_numdens(n, xc(1:n,1:3), radbdy(1:n),ndim, nbins, 1, n_mean(nsim,1:nbins), n_var(nsim,1:nbins), vol_bin(1:nbins))
    
    close(unit1, status='keep')

    OPEN(unit1, FILE= TRIM(RUN_NAME)//'_velpdf1_final.dat', status='replace')
    OPEN(unit2, FILE= TRIM(RUN_NAME)//'_velpdf2_final.dat', status='replace')
    if(ibidisperse) then
       do idim=1,ndim
          if(idim.eq.1) then 
             ftemp1(1:nbody1) = vx(1:nbody1)
             ftemp2(1:nbody2) = vx(nbody1+1:n)
          elseif(idim.eq.2) then 
             ftemp1(1:nbody1) = vy(1:nbody1)
             ftemp2(1:nbody2) = vy(nbody1+1:n)
          elseif(idim.eq.3) then 
             ftemp1(1:nbody1) = vz(1:nbody1)
             ftemp2(1:nbody2) = vz(nbody1+1:n)
          end if
          
          CALL histogram(ftemp1,wt1,nbody1,nhbins,fmin,fmax,hist(1:nhbins))
          CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,unit1,real(idim,prcn),1.d0)
          CALL histogram(ftemp2,wt2,nbody2,nhbins,fmin,fmax,hist(1:nhbins))
          CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,unit2,real(idim,prcn),1.d0)
       end do
    end if
    
    close(unit1, status='keep')
    close(unit2, status='keep')
    IF(nsim.LT.mis) goto 2001 

    rho_est_avg = zero 
    gofr_avg = zero 
    nmean_avg = zero
    nvar_avg = zero

    do j=1,nrbins
       rho_est_avg(j) = rho_est_avg(j) + 1./dble(nsim)*sum(rho_est(1:nsim,j))
       gofr_avg(j) = gofr_avg(j) + 1./dble(nsim)*sum(gofr(1:nsim,j))
    end do

    do j=1,nbins
       nmean_avg(j) = nmean_avg(j) + 1./dble(nsim)*sum(n_mean(1:nsim,j))
       nvar_avg(j) = nvar_avg(j) + 1./dble(nsim)*sum(n_var(1:nsim,j))
    end do
    do j=1,nbins
       IF(nsim.ge.2) THEN 
          conf1=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* &
               & sum((n_mean(1:nsim,j) - nmean_avg(j))**2))
          conf2=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* & 
               & sum((n_var(1:nsim,j)- nvar_avg(j))**2))
       end IF
       
       write(phisq_unit,'(10(E20.10,1x))')vol_bin(j),&
            &nmean_avg(j),nvar_avg(j),  conf1, conf2 
    end do
    close(phisq_unit, status="keep")
    do j=1,nrbins
       IF(nsim.ge.2) THEN 
          conf1=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* &
               & sum((rho_est(1:nsim,j) - rho_est_avg(j))**2))
          conf2=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* & 
               & sum((gofr(1:nsim,j)- gofr_avg(j))**2))
       end IF
       
       write(gofunit,'(10(E20.10,1x))')rad_bin(j),rad_bin(j)/r_0&
            &,rho_est_avg(j),gofr_avg(j),  conf1, conf2 
!&,rho_est_avg(j),gofr_avg(j),  conf1, conf2 
    end do
    close(gofunit, status = "keep") 

    

!!$    WRITE(*,'(//'' **** END OF DYNAMICS **** '')')
!!$
!!$    WRITE(*,'(/'' FINAL COLLIDING PAIR '',2I5)') I, J

    ! ** CHECK FOR PARTICLE OVERLAPS **

    ! CALL CHECK ( SIGMA, OVRLAP, E )

    !IF ( OVRLAP ) THEN

    !   WRITE(*,'('' PARTICLE OVERLAP IN FINAL CONFIGURATION '')')

    !ENDIF

    ! ** WRITE OUT CONFIGURATION **

    !CALL WRITCN ( CNFILE )

    !  ** WRITE OUT INTERESTING INFORMATION **

    PVNKT1 = ACW / REAL ( N ) / 3.0 / T / TEMP
    EN = E / REAL ( N )
    ENKT = EN / TEMP
    T = T * SQRT ( TEMP ) / SIGMA
    RATE = REAL ( NCOLL ) / T
    TBC  = REAL ( N ) / RATE / 2.0
    
    
    WRITE(1000,*) 'ZONE'
    !write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
    !     &    ' "UX" '
    DO i=1,n
       WRITE(1000,'(4(2x,f12.8))')rx(i), ry(i), rz(i), rp(i)
    enddo
    
    CLOSE(1000, status="keep")
    PRINT*,'PLACE1'
    
!!$    WRITE(*,'('' FINAL TIME IS          '',F15.8)') T
!!$    WRITE(*,'('' COLLISION RATE IS      '',F15.8)') RATE
!!$    WRITE(*,'('' MEAN COLLISION TIME    '',F15.8)') TBC
!!$    WRITE(*, '('' FINAL E/NKT IS         '',F15.8)') ENKT
!!$    WRITE(*,'('' PV/NKT - 1 IS          '',F15.8)') PVNKT1
    DEALLOCATE(rx,ry,rz, vx,vy, vz, wt,coltim&
         &,rp,rx_0,ry_0,rz_0, rxt,ryt,rzt )
    DEALLOCATE( partnr)
    if(ibidisperse)DEALLOCATE(wt1,wt2,ftemp1,ftemp2)
    PRINT*,'PLACE 2'
    
    RETURN 

    !STOP
  END SUBROUTINE PERFORM_HARD_SPHERE


  !***************************************************************
  SUBROUTINE UPLIST ( SIGMA, I )

    !  COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
    !  COMMON / BLOCK2 / COLTIM, PARTNR

    ! *******************************************************************
    ! ** LOOKS FOR COLLISIONS WITH ATOMS J > I                         **
    ! *******************************************************************

    double precision        timbig
    parameter ( timbig = 1.0e6 )

    integer, INTENT(IN) ::     i
    double precision, INTENT(IN)    ::    sigma

    integer     j
    double precision rxi, ryi, rzi, rxij, ryij, rzij
    double precision vxi, vyi, vzi, vxij, vyij, vzij
    double precision rijsq, vijsq, bij, discr, sigsq, tij

    double precision  rxti,ryti,rzti,rxtj,rytj,rztj,rijtsq
    double precision  rxtij, rytij, rztij, w

    double precision  rpijsq

    !     *******************************************************************

    if ( i .eq. n ) return

    coltim(i) = timbig
    rxi = rx(i)
    ryi = ry(i)
    rzi = rz(i)
    vxi = vx(i)
    vyi = vy(i)
    vzi = vz(i)

    do j = i + 1, n

       rxij = rxi - rx(j)
       ryij = ryi - ry(j)
       rzij = rzi - rz(j)
       rxij = rxij - anint ( rxij )
       ryij = ryij - anint ( ryij )
       rzij = rzij - anint ( rzij )
       vxij = vxi - vx(j)
       vyij = vyi - vy(j)
       vzij = vzi - vz(j)
       bij  = rxij * vxij + ryij * vyij + rzij * vzij

       if ( bij .lt. 0.d0 ) then

          rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2
          vijsq = vxij ** 2 + vyij ** 2 + vzij ** 2

          rpijsq  = (rp(i)+rp(j))**2

          discr = bij ** 2 - vijsq * ( rijsq - rpijsq )

          if ( discr .gt. 0.d0 ) then

             tij = ( -bij - sqrt ( discr ) ) / vijsq

             !write(*,*)coltemp,i, 'coltim...j..',coltim(j),j

             if ( tij .lt. coltim(i)) then
                if(tij.le.0.d0)then
                   write(*,*)'TIJ < 0 !!!!!!in uplist ..'
                   !write(*,*)vx(i),vy(i),vz(i)
                   !write(*,*)vx(j),vy(j),vz(j)
                   !write(*,*)rx(i),ry(i),rz(i)
                   !write(*,*)rx(j),ry(j),rz(j)
                   !write(*,*)rp(i),rp(j)
                   write(*,*)'coll.',coll,t!, tij, i, j, bij
                   write(*,*)rpijsq,rijsq, sqrt(rijsq)-sqrt(rpijsq)
                   overlap_uplist = .true. 
                   !call bump(sigma, i,j,w)
                   !bump_in_uplist = .TRUE.
                   partnr(i) = j
                   goto 40
                endif
                coltim(i) = tij
                partnr(i) = j
             endif
40           continue
          endif
       endif
    end do
1000 continue 
    return
  END SUBROUTINE UPLIST


  ! *******************************************************************
  SUBROUTINE BUMP ( SIGMA, I, J, W )


    ! *******************************************************************
    ! ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
    ! **                                                               **
    ! ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
    ! ** THE ROUTINE ALSO COMPUTES COLLISIONAL VIRIAL W.               **
    ! *******************************************************************


    integer     i, j
    double precision        sigma, w, m1, m2

    double precision        rxij, ryij, rzij, factor
    double precision  side, vxij,vyij,vzij, delvx, delvy, delvz, sigsq
    double precision  vijsq, rest0

    !     *******************************************************************

    sigsq = sigma ** 2

    m1 = 4./3.*pi*rp(i)**3.*rhod
    m2 = 4./3.*pi*rp(j)**3.*rhod
!!$    IF(COLL.EQ.1) THEN 
!!$       write(*,*)'Masses....',m1,m2
!!$       write(*,*)'Restitution coefficient..',rest
!!$    end IF
    
    rxij = rx(i) - rx(j)
    ryij = ry(i) - ry(j)
    rzij = rz(i) - rz(j)
    rxij = rxij - anint ( rxij )
    ryij = ryij - anint ( ryij )
    rzij = rzij - anint ( rzij )

    vxij = vx(i)-vx(j)
    vyij = vy(i)-vy(j)
    vzij = vz(i)-vz(j)

    vijsq = vxij ** 2 + vyij ** 2 + vzij ** 2 

    factor = ( rxij * ( vx(i) - vx(j) ) + ryij * ( vy(i) - vy(j) ) &
         &+rzij * ( vz(i) - vz(j) ) ) / (rp(i)+rp(j))**2 !sigsq 
    IF(factor.GT.0.d0)  THEN 
       WRITE(*,*)'ACHTUNG'
       WRITE(*,*) 'factor in bump GT zero = ', factor
       !STOP
    ENDIF
    
    delvx = - factor * rxij
    delvy = - factor * ryij
    delvz = - factor * rzij

    vx(i) = vx(i) + delvx * m2 *(1.d0+ rest)/(m1+m2)
    vx(j) = vx(j) - delvx * m1 *(1.d0+ rest)/(m1+m2)

    vy(i) = vy(i) + delvy * m2 *(1.d0+ rest)/(m1+m2)
    vy(j) = vy(j) - delvy * m1 *(1.d0+ rest)/(m1+m2)

    vz(i) = vz(i) + delvz * m2 *(1.d0+ rest)/(m1+m2)
    vz(j) = vz(j) - delvz * m1 *(1.d0+ rest)/(m1+m2)

    w = delvx * rxij + delvy * ryij + delvz * rzij

    return
  END SUBROUTINE BUMP

  subroutine computeener(nsp1, nsp2)

    integer, INTENT(in):: nsp1,nsp2
    Integer :: i, npart

    Double precision :: m, mtot, temp

    npart = nsp1+nsp2

    ener = 0.d0 
    ener1 = 0.d0 
    ener2 = 0.d0 
    ener3 = 0.d0 
    ener_spec1 = 0.d0
    ener_spec2 = 0.d0
    mass_weighted_var = 0.d0
    mtot = 0.d0

    do i=1,npart
       !if(rp(i).eq.radp1)then 
       m = 4./3.*pi*rp(i)**3.*rhod 
       temp =  (vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)) 
       if(i.le.nsp1)then
          ener_spec1 = ener_spec1 + temp
       else
          ener_spec2 = ener_spec2 + temp
       end if
       mass_weighted_var =  mass_weighted_var + m*temp
       mtot = mtot + m

       ener = ener + temp
       ener1 = ener1 + vx(i)*vx(i) 
       ener2 = ener2 + vy(i)*vy(i) 
       ener3 = ener3 + vz(i)*vz(i) 
    enddo
    ener_spec1 = ener_spec1/nsp1
    if(nsp2.gt.0) ener_spec2 = ener_spec2/nsp2
    mass_weighted_var =  mass_weighted_var/mtot

    ener = ener*0.5d0/dble(npart) 
    ener1 = ener1/dble(npart) 
    ener2 = ener2/dble(npart) 
    ener3 = ener3/dble(npart) 

    return

  end subroutine computeener




  SUBROUTINE DNLIST ( SIGMA, J )

    !  COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
    !  COMMON / BLOCK2 / COLTIM, PARTNR

    ! *******************************************************************
    ! ** LOOKS FOR COLLISIONS WITH ATOMS I < J                         **
    ! *******************************************************************

    INTEGER     N
    PARAMETER ( N = 108 )

    REAL(prcn)        TIMBIG
    PARAMETER ( TIMBIG = 1.E10 )

    INTEGER     J
    REAL(prcn)        SIGMA
    REAL(prcn)        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
    REAL(prcn)        COLTIM(N)
    INTEGER     PARTNR(N)

    INTEGER     I
    REAL(prcn)        RXJ, RYJ, RZJ, RXIJ, RYIJ, RZIJ
    REAL(prcn)        VXJ, VYJ, VZJ, VXIJ, VYIJ, VZIJ
    REAL(prcn)        RIJSQ, VIJSQ, BIJ, TIJ, DISCR, SIGSQ

    ! *******************************************************************

    IF ( J .EQ. 1 ) RETURN

    SIGSQ = SIGMA ** 2
    RXJ = RX(J)
    RYJ = RY(J)
    RZJ = RZ(J)
    VXJ = VX(J)
    VYJ = VY(J)
    VZJ = VZ(J)

    DO I = 1, J - 1

       RXIJ = RX(I) - RXJ
       RYIJ = RY(I) - RYJ
       RZIJ = RZ(I) - RZJ
       RXIJ = RXIJ - ANINT ( RXIJ )
       RYIJ = RYIJ - ANINT ( RYIJ )
       RZIJ = RZIJ - ANINT ( RZIJ )
       VXIJ = VX(I) - VXJ
       VYIJ = VY(I) - VYJ
       VZIJ = VZ(I) - VZJ
       BIJ  = RXIJ * VXIJ + RYIJ * VYIJ + RZIJ * VZIJ

       IF ( BIJ .LT. 0.0 ) THEN
          RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2
          VIJSQ = VXIJ ** 2 + VYIJ ** 2 + VZIJ ** 2
          DISCR = BIJ ** 2 - VIJSQ * ( RIJSQ - SIGSQ )
          IF ( DISCR .GT. 0.0 ) THEN
             TIJ = ( - BIJ - SQRT ( DISCR ) ) / VIJSQ
             IF ( TIJ .LT. COLTIM(I) ) THEN
                COLTIM(I) = TIJ
                PARTNR(I) = J
             ENDIF
          ENDIF
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE DNLIST

  

  SUBROUTINE CHECK ( SIGMA, OVRLAP, E )

    !  COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
    !  COMMON / BLOCK2 / COLTIM, PARTNR

    ! *******************************************************************
    ! ** TESTS FOR PAIR OVERLAPS AND CALCULATES KINETIC ENERGY.        **
    ! *******************************************************************

    INTEGER     N
    PARAMETER ( N = 108 )

    REAL(prcn)        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
    REAL(prcn)        COLTIM(N)
    INTEGER     PARTNR(N)

    REAL(prcn)        SIGMA, E
    LOGICAL     OVRLAP

    INTEGER     I, J
    REAL(prcn)        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ, RIJSQ, SIGSQ, RIJ
    REAL(prcn)        TOL
    PARAMETER ( TOL = 1.0E-4 )

    ! *******************************************************************

    SIGSQ  = SIGMA ** 2
    OVRLAP = .FALSE.
    E = 0.0

    DO I = 1, N - 1
       RXI = RX(I)
       RYI = RY(I)
       RZI = RZ(I)
       DO J = I + 1, N
          RXIJ = RXI - RX(J)
          RYIJ = RYI - RY(J)
          RZIJ = RZI - RZ(J)
          RXIJ = RXIJ - ANINT ( RXIJ )
          RYIJ = RYIJ - ANINT ( RYIJ )
          RZIJ = RZIJ - ANINT ( RZIJ )
          RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2
          IF ( RIJSQ .LT. SIGSQ ) THEN
             RIJ = SQRT ( RIJSQ / SIGSQ )
             WRITE(*,'('' I,J,RIJ/SIGMA = '',2I5,F15.8)') I, J, RIJ
             IF ( ( 1.0 - RIJ ) .GT. TOL ) OVRLAP = .TRUE.
          ENDIF
       ENDDO
    ENDDO

    DO I = 1, N
       E = E + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2
    ENDDO

    E = 0.5 * E

    RETURN
  END SUBROUTINE CHECK
end MODULE hard_sphere



    
