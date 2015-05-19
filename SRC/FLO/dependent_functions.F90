Module dependent_functions
#include "ibm.h"
	use precision
	use constants 
	use global_data
	!new use fftw_interface
	use fftw3_interface
	use parallel
	!use nlarrays
	implicit none
Contains
	subroutine set_interpolation_stencil(pc, ib,ie,jb,je,kb,ke, isch, ordernew)
		use global_data, Only : intx_per, inty_per, intz_per
		IMPLICIT NONE 
		integer, dimension(3), intent(in):: pc
		integer, intent(out):: ib, ie, jb, je, kb, ke 
		integer, OPTIONAL :: ordernew
		CHARACTER*5, intent(in) :: isch 
		integer :: ob2rtmp, ob2ltmp, ordertemp, ordernewtmp
		integer :: im, jm, km

		!new im = mxf
		im = local_ni(1)
		jm = local_ni(2)
		km = local_ni(3)

		SELECT CASE(isch)
		CASE('csi')
			ob2rtmp = ob2r
			ob2ltmp = ob2l
			ordertemp = order
			!!$if (order.NE.3.AND.(pc(1).EQ.1.OR.pc(1).EQ.im&
			!!$&-1.OR.pc(2).EQ.1.OR.pc(2).EQ.jm&
			!!$&-1.OR.pc(3).EQ.1.OR.pc(3).EQ.km-1)) order = 3
			!To make the order at boundary cells for csi equal to 3. 
			ob2l = (order+1)/2
			ob2r = order/2 
			!print*,'ob2l = ',ob2l,ob2r
		    
			if (.not.intx_per) then 
				ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
				ie = Min(im,pc(1) + ob2r)
				IF (ib.EQ.1 ) ie = ib + order - 1
				IF (ie.EQ.im) ib = ie - order + 1
			else 
				if (pc(1).lt.0) then
					ib = pc(1) - ob2r
					ie = pc(1) + (ob2l - 1) !periodic
				else
					ib = pc(1) - (ob2l - 1) !periodic
					ie = pc(1) + ob2r
				endif
			endif

			if (.not.inty_per) then
				jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
				je = Min(jm,pc(2) + ob2r)
				IF (jb.EQ.1 ) je = jb + order - 1
				IF (je.EQ.jm) jb = je - order + 1
			else
				if (pc(2).lt.0) then
					jb = pc(2) - ob2r
					je = pc(2) + (ob2l - 1) !periodic
				else
					jb = pc(2) - (ob2l - 1) !periodic
					je = pc(2) + ob2r
				endif
			endif
       
			if (.not.intz_per) then 
				kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
				ke = Min(km,pc(3) + ob2r)
				IF (kb.EQ.1 ) ke = kb + order - 1
				IF (ke.EQ.km) kb = ke - order + 1
			else
				if (pc(3).lt.0) then
					kb = pc(3) - ob2r
					ke = pc(3) + (ob2l - 1) !periodic
				else
					kb = pc(3) - (ob2l - 1) !periodic
					ke = pc(3) + ob2r
				endif
			endif

			ob2r =  ob2rtmp 
			ob2l = ob2ltmp
			ordernewtmp = order

			!print*,'ib,ie .... in processing =', ib,ie,jb,je,kb,ke,&
			!     & ordernewtmp
			!print*, 'pc = ',pc(1),pc(2),pc(3)!
			order = ordertemp !reset the order
		CASE('lpi')
			!print*, 'order in set stencil = ', order
			ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
			ie = Min(im,pc(1) + ob2r)
			if (.not.intx_per) then 
				IF (ib.EQ.1 ) ie = ib + order - 1
				IF (ie.EQ.im) ib = ie - order + 1
			else 
				IF (ib.EQ.1 ) ib = ie - order + 1
				IF (ie.EQ.im) ie = ib + order - 1
			endif

			jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
			je = Min(jm,pc(2) + ob2r)
			if (.not.inty_per) then
				IF (jb.EQ.1 ) je = jb + order - 1
				IF (je.EQ.jm) jb = je - order + 1
			else
				IF (jb.EQ.1 ) jb = je - order + 1
				IF (je.EQ.jm) je = jb + order - 1
			endif

			kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
			ke = Min(km,pc(3) + ob2r)
			if (.not.intz_per) then 
				IF (kb.EQ.1 ) ke = kb + order - 1
				IF (ke.EQ.km) kb = ke - order + 1
			else
				IF (kb.EQ.1 ) kb = ke - order + 1
				IF (ke.EQ.km) ke = kb + order - 1
			endif
			ordernewtmp = order
			end SELECT
		IF (PRESENT(ordernew)) ordernew = ordernewtmp
	end subroutine set_interpolation_stencil


	subroutine set_interpolation_pstencil(pc, ib,ie,jb,je,kb,ke, isch, ordernew)
		use global_data, Only :intx_per, inty_per,intz_per
		IMPLICIT NONE 
		integer, dimension(3), intent(in):: pc
		integer, intent(out):: ib, ie, jb, je, kb, ke 
		integer, OPTIONAL :: ordernew
		CHARACTER*5, intent(in) :: isch 
		integer :: im, jm, km
		integer :: ob2rtmp, ob2ltmp, ordertemp, ordernewtmp
		!new im = mxf - 1
		im = local_ni(1)
		jm = local_ni(2)
		km = local_ni(3)

		SELECT CASE(isch)
		CASE('csi')

			ob2rtmp = ob2r
			ob2ltmp = ob2l
			ordertemp = order
!!$       if (order.NE.3.AND.(pc(1).EQ.1.OR.pc(1).EQ.im&
!!$            &-1.OR.pc(2).EQ.1.OR.pc(2).EQ.jm&
!!$            &-1.OR.pc(3).EQ.1.OR.pc(3).EQ.km-1)) order = 3
			!To make the order at boundary cells for csi equal to 3. 
			ob2l = (order+1)/2
			ob2r = order/2 
			!print*,'ob2l = ',ob2l,ob2r

			if (.not.intx_per) then 
				ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
				ie = Min(im,pc(1) + ob2r)
				IF (ib.EQ.1 ) ie = ib + order - 1
				IF (ie.EQ.im) ib = ie - order + 1
			else 
				if (pc(1).lt.0) then
					ib = pc(1) - ob2r
					ie = pc(1) + (ob2l - 1) !periodic
				else
					ib = pc(1) - (ob2l - 1) !periodic
					ie = pc(1) + ob2r
				endif
			endif

			if (.not.inty_per) then
				jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
				je = Min(jm,pc(2) + ob2r)
				IF (jb.EQ.1 ) je = jb + order - 1
				IF (je.EQ.jm) jb = je - order + 1
			else
				if (pc(2).lt.0) then
					jb = pc(2) - ob2r
					je = pc(2) + (ob2l - 1) !periodic
				else
					jb = pc(2) - (ob2l - 1) !periodic
					je = pc(2) + ob2r
				endif
			endif
		    
		    
			if (.not.intz_per) then 
				kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
				ke = Min(km,pc(3) + ob2r)
				IF (kb.EQ.1 ) ke = kb + order - 1
				IF (ke.EQ.km) kb = ke - order + 1
			else
				if (pc(3).lt.0) then
					kb = pc(3) - ob2r
					ke = pc(3) + (ob2l - 1) !periodic
				else
					kb = pc(3) - (ob2l - 1) !periodic
					ke = pc(3) + ob2r
				endif
			endif

			ob2r =  ob2rtmp 
			ob2l = ob2ltmp
			ordernewtmp = order

			!print*,'ib,ie .... in processing =', ib,ie,jb,je,kb,ke,&
			!     & ordernewtmp
			!print*, 'pc = ',pc(1),pc(2),pc(3)!
			order = ordertemp !reset the order
		CASE('lpi')
			!print*, 'order in set stencil = ', order
			ib = MAX(1,  pc(1) - (ob2l - 1)) !non-periodic
			ie = Min(im, pc(1) + ob2r)
			if (.not.intx_per) then 
				IF (ib.EQ.1 ) ie = ib + order - 1
				IF (ie.EQ.im) ib = ie - order + 1
			else 
				IF (ib.EQ.1 ) ib = ie - order + 1
				IF (ie.EQ.im) ie = ib + order - 1
			endif
		    
			jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
			je = Min(jm,pc(2) + ob2r)
			if (.not.inty_per) then
				IF (jb.EQ.1 ) je = jb + order - 1
				IF (je.EQ.jm) jb = je - order + 1
			else
				IF (jb.EQ.1 ) jb = je - order + 1
				IF (je.EQ.jm) je = jb + order - 1
			endif
		    
			kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
			ke = Min(km,pc(3) + ob2r)
			if (.not.intz_per) then 
				IF (kb.EQ.1 ) ke = kb + order - 1
				IF (ke.EQ.km) kb = ke - order + 1
			else
				IF (kb.EQ.1 ) kb = ke - order + 1
				IF (ke.EQ.km) ke = kb + order - 1
			endif
			ordernewtmp = order
		end SELECT
		if (PRESENT(ordernew)) ordernew = ordernewtmp
	end subroutine set_interpolation_pstencil

#if 0
  subroutine calc_realdiv(divur)
    !use nlarrays, Only : uf1
    implicit none 

    integer :: i,j,k

    real(prcn), dimension(:,:,:), intent(out) ::  divur
    do i=1,mx1,1
       do j=1,my2,1
          do k=1,mz,1
             uf1(j,k)=(u(i+1,j,k,1)-u(i,j,k,1))/dx
             uf1(j,k)=uf1(j,k)+half*wy(j)*(u(i,j,k,2)+u(i+1,j,k,2))
             uf1(j,k)=uf1(j,k)+half*wz(k)*(u(i,j,k,3)+u(i+1,j,k,3))
          enddo
       enddo

       call ff2cr(uf1,divur(i,:,:))
    enddo
  end subroutine calc_realdiv
#endif

#if 0
  subroutine calc_forreal(fr)
    implicit none 

    real(prcn), dimension(:,:,:,:), intent(out) :: fr
    integer :: i,n
    do n=1,ndim
       do i=1,mxf

          call ff2cr(ff(i,:,:,n),fr(foffset+i,:,:,n))
          
       enddo			! loop over i
    enddo

  end subroutine calc_forreal


  subroutine calc_velreal_mxf(velr)
    implicit none 

    real(prcn), dimension(:,:,:,:), intent(out) :: velr
    integer :: i,j,k, n ,ii

    do n=1,ndim

       do i=1,mxf,1
          ii = i+foffset
          call ff2cr(u(ii,:,:,n),velr(i,:,:,n))
          velr(i,:,:,n) = velr(i,:,:,n) + umean(n)
       enddo
    enddo

  end subroutine calc_velreal_mxf


  subroutine calc_velreal(velr)
    implicit none 

    real(prcn), dimension(:,:,:,:), intent(out) :: velr
    integer :: i,j,k, n, idim


#if PARALLEL
	do idim = 1, ndim
		i = 1
		VECSendRECV(u(i,1,1,idim),1,ucslice,fromproc,1,u(nx+1,1,1,idim),1,toproc,1,comm_cart_2d,status)
		i = nx
		VECSendRECV(u(i,1,1,idim),1,ucslice,toproc,0,u(0,1,1,idim),1,fromproc,0,comm_cart_2d,status)
	enddo
#endif


    do n=1,ndim
#if PARALLEL
       do i=0,nx+1
#else
       do i=1,nx !mx,1       
#endif
          call ff2cr(u(i,:,:,n),velr(i,:,:,n))
          velr(i,:,:,n) = velr(i,:,:,n) + umean(n)
       enddo
#if !PARALLEL
			velr(mx,:,:,n) = velr(1,:,:,n)
#endif
    enddo

  end subroutine calc_velreal
#endif

	subroutine calc_pressure
		use nlmainarrays, only : pbcp
		implicit none 
		call fftwc2r(p, pbcp(1:local_ni(1), 1:local_ni(2), 1:local_ni(3)))
		call communicate_in_gohst_domain(pbcp)
	end subroutine calc_pressure


!  subroutine calc_realvort
!    use nlarrays , Only : uf1
!    implicit none 
!
!    integer :: i,j,k
!  end subroutine calc_realvort
  
  subroutine writepres_sph!(pr,ur) 
    use general_funcs
    use interpolation

    implicit none 
    integer :: unitno, m, l
    integer, SAVE :: count_routine = 0 
    real(prcn) :: stagangle 
    CHARACTER*100 :: FILENAME1
    logical :: filexist, isopen 
    count_routine  = count_routine + 1
	
    if (count_routine.eq.3) count_routine = 1

 3020 FORMAT(I1)
    write (FILENAME1, 3020) count_routine	
    FILENAME1 = TRIM(RUN_NAME)//'_cpvstheta_'//TRIM(FILENAME1)//'.dat'
   inQUIRE(FILE=FILENAME1,EXIST=filexist,OPENED=isopen)
       
    unitno = getnewunit(minunitno, maxunitno)
    IF (.NOT.filexist) then
       
    OPEN(unitno,FILE=FILENAME1, status='new')
      elseif (filexist.AND..NOT.isopen) then
    open(unit=unitno,file=FILENAME1,status='replace')
    endif

    write (unitno, '("# generated at time = ",g12.5)') t
    do m=1,nbody              ! loop over bodies
       write (unitno,*) 'zone'
       do l=1,nbnd
          if (xs(3,l).eq.zero)then
             if (xs(2,l).ge.zero)then
                if (xs(1,l).lt.zero)then
                   stagangle = abs(atan(xs(2,l)/xs(1,l)))
                elseif (xs(1,l).gt.zero)then
                   stagangle = atan(xs(2,l)/xs(1,l))
                   stagangle = pi-stagangle
                elseif (xs(1,l).eq.zero)then
                   stagangle = pi/2.
                endif
                
!                write (unitno,101)stagangle*180./pi,coeffp(m,l)
                !     end loop over all bnd points
             endif
          endif
       enddo
101    FORMAT(10(2x,e20.12))
       
       !     end loop over all bodies
    enddo
       close(unitno, status = "keep")
  end subroutine writepres_sph
  
    
  subroutine deterpressph!(pr,ur) 
    use general_funcs
    use interpolation

    use nlmainarrays, Only : ur=>ubcp, nlr=>nlbcp, onlr=>onlbcp,pr=> pbcp
    use bcsetarrays, Only : ppr
    !use boundary_condition
    implicit none 

    !real(prcn), intent(in) , dimension(:,:,:) ::  pr
    !real(prcn), intent(in) , dimension(:,:,:,:) ::  ur


    real(prcn) ::  xl(ndim), xpb(ndim)
    real(prcn) ::  ul(ndim),unorm(ndim),utang(ndim), ppll(ndim)
    real(prcn) ::  snorm(ndim), stang(ndim)
    real(prcn) ::  unmag,tx,ty,tz, dpdthetamag

    integer :: m,l,n,i,j,k,d,pcell(3), ii,jj,kk

    real(prcn) ::  rad,pl, dfll(ndim), nll(ndim), onll(ndim)
    real(prcn) ::  stagangle
    integer :: isp
    integer :: is(ndim)
    integer :: ib,ie, jb,je, kb,ke , onew
    open(unit=45,file='presphperip.dat',status='unknown')


   print*,'in DETER PRESS: NBODY = ', nbody, nbnd 
    write (100,*)'Zone'

    do m=1,nbody              ! loop over bodies
       write (45,*) 'zone'
       do l=1,nbnd

          if (xs(3,l).eq.zero)then
             if (xs(2,l).ge.zero)then
                rad=0.0
                do n=1,ndim

                   xl(n)=xc(m,n)+xs(n,l)*radbdy(m)
                   is(n)=int(xl(n))

                   ul(n)=zero
                   rad=rad+(xs(n,l)*radbdy(m))**2.0
                enddo

                


                rad=dsqrt(rad)
                
                if (xs(1,l).lt.zero)then
                   stagangle = abs(atan(xs(2,l)/xs(1,l)))
                elseif (xs(1,l).gt.zero)then
                   stagangle = atan(xs(2,l)/xs(1,l))
                   stagangle = pi-stagangle
                elseif (xs(1,l).eq.zero)then
                   stagangle = pi/2.
                endif
                do n=1,ndim
                   snorm(n)=(xs(n,l)*radbdy(m))/rad
                enddo
                stang(1) = snorm(2)
                stang(2) = -snorm(1)
                stang(3) = zero

                write (100,'(4(F15.10,1x))')xs(1,l),xs(2,l),xs(3,l)&
                     &,stagangle

                pl = zero
                ppll = zero
                isp=inT(xl(1)-0.5)
                pcell(1) = isp 
                pcell(2:3) = is(2:3)
                xpb(1) = xl(1)-0.5
                xpb(2:3)=xl(2:3)
!                call interpolate_pdata(pcell,xpb, ppr,ppll,pl)
!!$                call set_interpolation_pstencil(pcell,ib,ie,jb,je,kb,ke&
!!$                     &,interp_scheme, onew) 
!!$                do k = 1, onew
!!$                   do j = 1, onew
!!$                      do i = 1, onew
!!$                         ii = ib+i-1
!!$                         jj = jb+j-1
!!$                         kk = kb+k-1
!!$                         gstencil(i,j,k,1) = ib+(i-1)
!!$                         gstencil(i,j,k,2) = jb+(j-1)
!!$                         gstencil(i,j,k,3) = kb+(k-1)
!!$	     
!!$                         if (ii.lt.1.and.intx_per) then
!!$                            !print*,'ii LT 1'
!!$                            ii = mxf+ii-1
!!$                         endif
!!$             
!!$                         if (ii.gt.mxf-1.and.intx_per) then 
!!$                            ii = ii-mxf +1
!!$                            !print*,'ii GT mXF-1'
!!$                            !READ(*,*)
!!$                         endif
!!$                         
!!$                         
!!$                         if (jj.lt.1) jj = my+jj
!!$                         if (jj.gt.my) jj = jj-my
!!$                         if (kk.lt.1) kk = mz+kk
!!$                         if (kk.gt.mz) kk = kk-mz
!!$                         prsten(i,j,k) = pr(ii,jj,kk)
!!$                      enddo
!!$                   enddo
!!$                enddo
!!$                call interpolator(gstencil(1:onew,1:onew,1:onew,1:ndim)&
!!$                     &,prsten(1:onew,1:onew,1:onew)&
!!$                     &,xpb(1:ndim),pl,onew, interp_scheme,weightp)
                dpdthetamag = zero
                do n = 1, ndim
                   dpdthetamag = dpdthetamag + stang(n)*ppll(ndim)
                enddo
                
!!$                write (45,101)stagangle*180./pi,pl/0.5/(upi(1)**2.)
                write (45,101)stagangle*180./pi,dpdthetamag/0.5/(uchar(1)**2.)
             endif
          endif
       enddo

       !     end loop over all bodies
    enddo

101 FORMAT(10(2x,e20.12))
    close(45)


  end subroutine deterpressph


	subroutine interpolate_udata(pc,pos,ib,ie,jb,je,kb,ke,ul,nll,onll,dfll,flag, ibody,l,onew)
		use general_funcs
		use interpolation
		use bcsetarrays, ONLY :  diffn 
		use nlmainarrays, Only : ur=>ubcp, nlr=>nlbcp, onlr=>onlbcp!,pr=> pbc 
		implicit none 
		integer, intent(in) :: pc(3)
		integer, intent(in) :: flag, ibody,l
		real(prcn), dimension(:), intent(in) :: pos
		integer, intent(out) :: ib, ie, jb,je,kb,ke, onew
		real(prcn), intent(out), dimension(:) :: ul,nll,onll,dfll !, ppll
		integer :: i, j,k, ii,jj,kk, n
		logical :: detect_error
		!!$real(prcn), dimension(:,:,:,:), intent(in) ::  ppr

		call set_interpolation_stencil(pc,ib,ie,jb,je,kb,ke,interp_scheme, onew) 
		!new if ((ib.lt.1.or.ie.gt.mxf) .and. .not.xperiodic) print*,'Error in i ....',ib,ie,pc,pos
		!if ((ib.lt.1.or.ie.gt.mx) .and. .not.xperiodic) print*,'Error in i ....',ib,ie,pc,pos
		detect_error = .false.
		do k = 1, onew
			do j = 1, onew
				do i = 1, onew
					ii = ib+i-1
					jj = jb+j-1
					kk = kb+k-1
					gstencil(i,j,k,1) = ib+(i-1)
					gstencil(i,j,k,2) = jb+(j-1)   
					gstencil(i,j,k,3) = kb+(k-1)

					if (debug_check) then
						if (flag==0) then
							if (ii<-1 .or. ii>local_ni(1)+2) detect_error = .true.
							if (jj<-1 .or. jj>local_ni(2)+2) detect_error = .true.
							if (kk<-1 .or. kk>local_ni(3)+2) detect_error = .true.
						elseif (flag==1) then
							if (ii<0 .or. ii>local_ni(1)+1) detect_error = .true.
							if (jj<0 .or. jj>local_ni(2)+1) detect_error = .true.
							if (kk<0 .or. kk>local_ni(3)+1) detect_error = .true.
						endif

						if (detect_error) then
							write (*,"(1A,1i6,1a,3d15.7,1a,3i6)") "ERROR DETECTED IN VELOCITY POINT IN ID: ", myid, ", POSITION:", pos(:), &
																			& " ,INDEX: ", ii, jj, kk
							PARALLEL_FINISH()
							stop
						endif
					endif

					vsten(i,j,k,1:ndim) = ur(ii,jj,kk,1:ndim)
					if (flag.eq.1) then 
						nlsten(i,j,k,1:ndim) = nlr(ii,jj,kk,1:ndim)
						onlsten(i,j,k,1:ndim) = onlr(ii,jj,kk,1:ndim)
						dfsten(i,j,k,1:ndim) = diffn(ii,jj,kk,1:ndim)
					endif
				enddo
			enddo
		enddo

		call interpolator(gstencil(1:onew,1:onew,1:onew,1:3), vsten(1:onew,1:onew,1:onew,1:ndim), pos(1:ndim),ul(1:ndim),onew, interp_scheme,weightp) 

		if (flag.eq.1) then 
			do n = 1, ndim 
				nll(n) =  array_dot_product(nlsten(1:onew,1:onew,1:onew,n),weightp(1:onew,1:onew,1:onew)) 
				onll(n)=  array_dot_product(onlsten(1:onew,1:onew,1:onew,n),weightp(1:onew,1:onew,1:onew))
				dfll(n)=  array_dot_product(dfsten(1:onew,1:onew,1:onew,n),weightp(1:onew,1:onew,1:onew))
				!!$ppll(n)=  array_dot_product(ppgrsten(1:onew,1:onew,1:onew,n),weightp(1:onew,1:onew,1:onew))
			enddo
		endif
	end subroutine interpolate_udata

	subroutine interpolate_pdata(pc,pos,ppl,pl,l)
		use general_funcs
		use interpolation
		use nlmainarrays, Only : pr=> pbcp 
		use bcsetarrays, Only : ppr
		implicit none 
		integer, intent(in) :: pc(3),l
		!real(prcn), dimension(:,:,:,:), intent(in) ::  ppr
		!real(prcn), dimension(:,:,:), intent(in) ::  pr
		real(prcn), dimension(:), intent(in) :: pos
		real(prcn), intent(out), dimension(:) :: ppl
		real(prcn), intent(out) :: pl
		integer :: i, j,k, onew, ii,jj,kk, n,ib,ie,jb,je,kb,ke, send_error, rec_error
		logical :: detect_error

		call set_interpolation_pstencil(pc,ib,ie,jb,je,kb,ke,interp_scheme, onew) 
		detect_error = .false.
		do k = 1, onew
			do j = 1, onew
				do i = 1, onew
					ii = ib+i-1
					jj = jb+j-1
					kk = kb+k-1
					gstencil(i,j,k,1) = ib+(i-1)
					gstencil(i,j,k,2) = jb+(j-1)
					gstencil(i,j,k,3) = kb+(k-1)

					if (debug_check) then
						if (ii<0 .or. ii>local_ni(1)+1) detect_error = .true.
						if (jj<0 .or. jj>local_ni(2)+1) detect_error = .true.
						if (kk<0 .or. kk>local_ni(3)+1) detect_error = .true.

						if (detect_error) then
							write (*,"(1A,1i6,1a,3d15.7,1a,6i6)") "ERROR DETECTED IN PRESSURE POINT IN ID: ", myid, ", POSITION:", pos(:), &
																			& " ,INDEX: ", ib,ie,jb,je,kb,ke
							PARALLEL_FINISH()
							stop
						endif
					endif

					ppgrsten(i,j,k,1:ndim) = ppr(ii,jj,kk,1:ndim)
					prsten(i,j,k) = pr(ii,jj,kk)!-Dreal(p(ii,1,1))
				enddo
			enddo
		enddo
		call interpolator(gstencil(1:onew,1:onew,1:onew,1:ndim),ppgrsten(1:onew,1:onew,1:onew,1:ndim)&
			&,pos(1:ndim), ppl(1:ndim), onew, interp_scheme,weightp)
		pl = array_dot_product(prsten(1:onew,1:onew,1:onew), weightp(1:onew,1:onew,1:onew)) 
  end subroutine interpolate_pdata


	subroutine grid_nodes_insphere
		Implicit None 

		integer :: i,j,k, ii, idim, cor_min(ndim), cor_max(ndim), m, imin, imax, jmin, jmax, kmin, kmax, iphs
		real(prcn) :: xlr(ndim), xll(ndim), dist, volfracg(nphases)
		logical :: already_accessed
		integer :: numpart_accessed, part_neigh, gn_pneigh_index
		integer :: index_in(ndim), index_out(ndim)
		real(prcn) :: position_in(ndim), position_out(ndim)
		logical :: i_have_the_point

		integer(8) :: count_solid_tmp, count_fluid_tmp

		if (I_AM_NODE_ZERO) write (*,'(A)')'in GRID NODES in SPHERE'
		do iphs = 1, nphases
			phase_array(iphs)%volfracg = zero
			volfracg(iphs) = zero
		enddo

		if (.not.allocated(myid_particles)) allocate(myid_particles(nbody))
		myid_particles(:) = .false.


		!new do k = 1, mz
		!new 	do j = 1, my
		!new 		do i = 0, nx+1
		fluid_atijk(:,:,:) = .true.
		do k=1, local_ni(3)
			do j=1, local_ni(2)
				do i=1, local_ni(1)
					gnacc_part(i,j,k,1) = 0
					do ii = 2, maxgn_partneigh
						gnacc_part(i,j,k,II) = -1
					enddo
				enddo
			enddo
		enddo
		count_solid  = 0
		count_solid_tmp = 0
		count_fluid_tmp  = local_ni(1)*local_ni(2)*local_ni(3)

		do m=1, nbody
			do idim=1, ndim 
				!if (idim.eq.1) then 
				!	xlr(idim) = xc(m,idim)  + radbdy(m) !+ foffset
				!	xll(idim) = xc(m,idim)  - radbdy(m) !+ foffset
				!else 
					xlr(idim) = xc(m,idim)  + radbdy(m)
					xll(idim) = xc(m,idim)  - radbdy(m) 
				!endif
			enddo
    
			do idim = 1, ndim 
				cor_min(idim) = floor(xll(idim)) !ceiling(xll(idim))
				cor_max(idim) = ceiling(xlr(idim))!floor(xlr(idim)) 
			enddo
    
			imin = cor_min(1) - 2*radbdy(m)
			imax = cor_max(1) + 2*radbdy(m)
			jmin = cor_min(2) - 2*radbdy(m)
			jmax = cor_max(2) + 2*radbdy(m)
			kmin = cor_min(3) - 2*radbdy(m)
			kmax = cor_max(3) + 2*radbdy(m)

			do i=imin, imax
				do j=jmin, jmax
					do k=kmin, kmax
						index_in(1) = i
						index_in(2) = j
						index_in(3) = k

						position_in(1) = i
						position_in(2) = j
						position_in(3) = k

						i_have_the_point = point_in_this_domain(index_in,index_out,position_in,position_out)
						if (i_have_the_point) then

							myid_particles(m) = .true.
							!THE ABOVE LINE AND THE BOTTOM LINE ARE ADDED TO FIND OUT WHAT PARTICLES ARE IN PROXIMITY OR INSIDE A DOMAIN
							!TO SAVE TIME BY NOT LOOPING OVER ALL PARTICLES IN DRAG AND IB CALCULATIONS
							if (cor_min(1)<=i.and.i<=cor_max(1) .and. cor_min(2)<=j.and.j<=cor_max(2) .and. cor_min(3)<=k.and.k<=cor_max(3)) then
								numpart_accessed = gnacc_part(index_out(1),index_out(2),index_out(3),1)
								already_accessed = .false.
								do part_neigh = 2, numpart_accessed+1
									if (m.eq.gnacc_part(index_out(1),index_out(2),index_out(3),part_neigh)) already_accessed = .true.
								enddo
				          
								if (.not.already_accessed) then
									gnacc_part(index_out(1),index_out(2),index_out(3),1) = gnacc_part(index_out(1),index_out(2),index_out(3),1) + 1
									gn_pneigh_index = gnacc_part(index_out(1),index_out(2),index_out(3),1) + 1
									if (gn_pneigh_index.gt.maxgn_partneigh+1) then
										write (*,*) 'gn_pneigh_index =',gn_pneigh_index ,' > maxgn_partneigh =', maxgn_partneigh, &
													& ' FOR GRID NODE', index_out(1), index_out(2), index_out(3)
										PARALLEL_FINISH()
										stop
									endif
									gnacc_part(index_out(1),index_out(2),index_out(3),gn_pneigh_index) = m
								endif

								dist = (i - xc(m,1))**two + (j - xc(m,2))**two + (k - xc(m,3))**two 
								dist = dsqrt(dist)
								if ((dist - radbdy(m)).LE.small_number) then 
									if (fluid_atijk(index_out(1),index_out(2),index_out(3)).eq..false.) then 
										write (*,"(1a,5i8)") 'fluid_atijk ALREADY false AT I,J,K=', index_out(1),index_out(2),index_out(3), myid, m
									else
										fluid_atijk(index_out(1),index_out(2),index_out(3))  = .false.
										iphs = part_array(m)%iphs
										volfracg(iphs) = volfracg(iphs) + one
										count_solid_tmp = count_solid_tmp + 1
										count_fluid_tmp = count_fluid_tmp - 1
									endif
								endif
							endif
						endif
					enddo
				enddo
			enddo
		enddo
		call communicate_fluid_atijk

		GLOBAL_LONGINT_SUM(count_solid_tmp, count_solid, 1, comm_cart_2d)
		GLOBAL_LONGINT_SUM(count_fluid_tmp, count_fluid, 1, comm_cart_2d)
		count_total = count_solid + count_fluid

		do iphs = 1, nphases
			GLOBAL_DOUBLE_SUM(volfracg(iphs),phase_array(iphs)%volfracg,1,comm_cart_2d)
		enddo

		maxvolfrac = real(count_solid,prcn)/real(count_total,prcn)
		do iphs = 1, nphases
			phase_array(iphs)%volfracg = phase_array(iphs)%volfracg/real(count_total,prcn)
			if (I_AM_NODE_ZERO) write (*,"(1a,1i,1a,1d15.7)") 'PHASE ', iphs, ' VOLUME FRACTION = ', (phase_array(iphs)%volfracg)
		enddo
	end subroutine grid_nodes_insphere


	subroutine compute_mean_fluid_velocity
		IMPLICIT NONE
		integer :: m, iphs, n, i, j, k, ii, jj, kk, idim
		integer :: l, vcellb(ndim), ci, cj, ck, numpart_accessed, PNEIGH
		integer :: ib, ie, jb, je, kb, ke, onew
		real(prcn) ::  xl(ndim), ul(ndim), nll(ndim), dfll(ndim), onll(ndim), dist, signed_dist, total_signed_dist
		real(prcn) :: total_solid_volfrac, cell_volfraction
		real(prcn) :: ldom(ndim), tempr(ndim), total_fluid_vol, fluid_meanvel(ndim)

		total_solid_volfrac = zero
		total_fluid_vol = zero
		fluid_meanvel(:) = zero

		ldom(1) = real(global_n(1), prcn)
		ldom(2) = real(global_n(2), prcn)
		ldom(3) = real(global_n(3), prcn)

		do ck=1, local_ni(3)
			do cj=1, local_ni(2)
				do ci=1, local_ni(1)
					xl(1) = (ci-half)
					xl(2) = (cj-half)
					xl(3) = (ck-half)

					do n = 1, ndim
						vcellb(n) = floor(xl(n))
					enddo
             
					call interpolate_udata(vcellb,xl,ib,ie,jb,je,kb,ke,ul,nll,onll,dfll, 0,-1, l, onew) 

					cell_volfraction = zero
					total_signed_dist = zero

					do k = 1, onew 
						do j = 1, onew
							do i = 1, onew
								ii = ib+i-1
								jj = jb+j-1
								kk = kb+k-1
#if !PARALLEL
								!new if (ii.lt.1) ii = mxf+ii-1
								!new if (ii.gt.mxf-1) ii = ii-(mxf-1)
								!if (ii.lt.1) ii = mx+ii
								!if (ii.gt.mx) ii = ii-mx
#endif
								!if (jj.lt.1) jj = my+jj
								!if (jj.gt.my) jj = jj-my
								!if (kk.lt.1) kk = mz+kk
								!!if (kk.gt.mz) kk = kk-mz 
								!LOCAL_inDEX(ii)                

								numpart_accessed = gnacc_part(ii,jj,kk,1)
								!!$ write (*,*)'I AM HERE ::', ii, jj, kk, numpart_accessed
								signed_dist = large_number

								do PNEIGH = 2, numpart_accessed+1
									m = gnacc_part(ii,jj,kk,PNEIGH)
									TEMPR(1) = ii - XC(m,1)
									TEMPR(2) = jj - XC(m,2)
									TEMPR(3) = kk - XC(m,3)

									dist = zero
									do idim = 1, ndim
										if ((ABS(tempr(idim))).gt.ldom(idim)/two) then
											if (tempr(idim).lt.zero) then
												tempr(idim) = tempr(idim) + Ldom(idim)
											else
												tempr(idim) = tempr(idim) - Ldom(idim)
											endif
										endif
										dist = dist + tempr(idim)**2.d0
									enddo
									dist = dsqrt(dist)
									dist = dist - RADBDY(m)
									signed_dist = Min(dist, signed_dist)
								enddo
								if (.not.fluid_atijk(ii,jj,kk)) cell_volfraction = cell_volfraction - signed_dist
								total_signed_dist = total_signed_dist + ABS(signed_dist)
							enddo
						enddo
					enddo
					cell_volfraction = cell_volfraction/total_signed_dist
					fluid_meanvel(:) = fluid_meanvel(:) + ul(:)*(one-cell_volfraction)*dx**3.d0
					total_fluid_vol = total_fluid_vol + (one-cell_volfraction)*dx**3.d0
					total_solid_volfrac = total_solid_volfrac + cell_volfraction
				enddo
			enddo
		enddo
    
		total_solid_volfrac = total_solid_volfrac/real(global_n(1)*global_n(2)*global_n(3),prcn)
		fluid_meanvel(:)  = fluid_meanvel(:)/total_fluid_vol

		ufmean(:) = fluid_meanvel(:)
		!write (*,*)'TOTAL SOLID VOLUME FRACTION: ', total_solid_volfrac
		!write (*,*)'MEAN FLUID VEL NEW METHOD: ', fluid_meanvel(:)
	end subroutine compute_mean_fluid_velocity

	subroutine update_nrpr_array(m)
		IMPLICIT NONE
		integer, intent(in) :: m
		integer :: l,n, iphs
		real(prcn) :: xl(ndim), xli(ndim), xlo(ndim)


		iphs = 1!part_array(m)%iphs

		part_array(m)%nrpr_active = phase_array(iphs)%nrpr
		part_array(m)%drag_active = phase_array(iphs)%nrpr

		NULLIFY(bndarray)
		bndarray => phase_array(iphs)%bndpts
		nrpr = part_array(m)%nrpr_active

!if (m==2) debug_check = .true.

		do l=1,nrpr
			do n=1,ndim
				xl(n)  = xc(m,n) + bndarray(n,l) * radbdy(m)
				xli(n) = xc(m,n) + bndarray(n,l) * radibdy(m)
				xlo(n) = xc(m,n) + bndarray(n,l) * radobdy(m)
			enddo

!if (debug_check) then
!	write (*,"(3d15.7)") radibdy(m), radbdy(m), radobdy(m)
!endif

			call check_external_pt(xlo, xl, m, l)

			if (.not.(part_array(m)%if_rev(l)))  part_array(m)%nrpr_active = part_array(m)%nrpr_active - 1
			if (.not.(part_array(m)%if_drag(l))) part_array(m)%drag_active = part_array(m)%drag_active - 1


!if (debug_check) then
!	write (*,*) part_array(m)%nrpr_active, part_array(m)%drag_active
!endif


		enddo
  end subroutine update_nrpr_array

	subroutine check_external_pt(xloin, xlin, m, l)
		use dem_mod
		implicit none 
		real(prcn), dimension(:), intent(in) :: xloin, xlin
		integer, intent(in) :: m, l
		!!$ logical, intent(OUT) :: REVERSE
		!new logical :: PER
		real(prcn), dimension(ndim) :: xlo, xl
		integer :: ip, PNO, IDIM

!integer :: i,j

		real(prcn) :: disto(ndim), dist(ndim), rado, rad, rmax(ndim)
		!!$ REVERSE = .true.

		part_array(m)%if_rev(l)  = .true.
		part_array(m)%if_drag(l) = .true.
		!new PER = .false.
		xlo = xloin
		xl  = xlin

		rmax(:) = dble(global_n(:))/2

		!new if (xlo(1).lt.one) xlo(1) = mxf+xlo(1)
		!new if (xlo(1).gt.mxf) xlo(1) = xlo(1)-mxf
		!new if (xlo(2).lt.one) then 
		!new	xlo(2) = my+one+xlo(2)
		!new	!print*, 'xlo2 LT zero'
		!new	PER = .true.
		!new endif
		!new if (xlo(2).gt.my) then 
		!new	xlo(2) = xlo(2)-my-one
		!new	!print*, 'xlo2 GT MY'
		!new	PER = .true.
		!new endif
		!newif (xlo(3).lt.one) then 
		!new 	xlo(3) = mz+one+xlo(3)
		!new 	PER = .true.
		!new endif
		!new if (xlo(3).gt.mz) then 
		!new 	xlo(3) = xlo(3)-mz-one
		!new 	PER = .true.
		!new endif
    
#if 0
do i=1, nbody
	write (*,*) "MY id = , neighbours", i, neighbours(i,1)+1
	write (*,*) neighbours(i,:)
enddo
read (*,*)
#endif
	


		if (neighbours(m,1).GT.1) then 
			do ip = 2, neighbours(m,1)+1
				PNO = neighbours(m,ip)
				!print*,'PNIO = ', PNO, neighbours(m,neighbours(m,1)+1)
				if (pno/=m) then
					do idim = 1, ndim
						disto(idim) = abs(xc(pno,idim)-xlo(idim))
						if (disto(idim)>rmax(idim)) disto(idim) = 2*rmax(idim) - disto(idim)

						dist(idim) = abs(xc(pno,idim)-xl(idim))
						if (dist(idim)>rmax(idim)) dist(idim) = 2*rmax(idim) - dist(idim)
					enddo
					rado = sqrt(dot_product(disto(:),disto(:)))
					rad  = sqrt(dot_product(dist(:),dist(:)))

					if (rado < radbdy(pno)) then
						!write (*,*) 1
						part_array(m)%if_rev(l) = .false.
					endif

					if (rad < radbdy(pno)) then
						!write (*,*) 2
						part_array(m)%if_drag(l) = .false.
					endif

#if 0
if (debug_check) then
	write (*,*) "part, neighb = ", m, pno
	write (*,"(1a,3d15.7)") "xlo   = ", xlo
	write (*,"(1a,3d15.7)")  "xl    = ", xl

	write (*,"(1a,3d15.7)")  "disto = ", disto
	write (*,"(1a,3d15.7)")  "dist  = ", dist

	write (*,"(1a,3d)")  "rado  = ", rado, radbdy(pno)
	write (*,"(1a,3d)")  "rad   = ", rad, radbdy(pno)

	write (*,*) part_array(m)%if_rev(l), part_array(m)%if_drag(l)
	read(*,*)
endif
#endif
					if ((.not.part_array(m)%if_rev(l)) .and. (.not.part_array(m)%if_drag(l))) exit
				endif
			enddo
		endif
	end subroutine check_external_pt



	subroutine compute_new_timestep(rks)
		use mypost_process, only : reynolds_stress_tensor, compute_sijsij
		use init_turb, only : calc_velreal
		use nlmainarrays, Only : ubcp, pbcp
		use dem_mod, only : is_mobile
		IMPLICIT NONE
		integer, intent(in) :: rks
		real(prcn) ::  u_max, v_max, w_max, umax_tmp, udiff
		real(prcn) ::  umax_loc, vmax_loc, wmax_loc  
		real(prcn) :: mixmeanslip(ndim),umean_temp(ndim), usmean_temp(ndim), ufmean_temp(ndim), &
						& umean_temploc(ndim), usmean_temploc(ndim), ufmean_temploc(ndim), mesh_veltemp(ndim)
    
		integer :: idim, i, j, k, iphs, partstart, partend, m
		real(prcn) :: dt_ratio, max_part_vel, part_vel_tmp(ndim), phase_mass(nphases)

		real(prcn) cpu0, cpu1

		if (I_AM_NODE_ZERO) call CPU_TIME (CPU0)


		dtnew = large_number
		u_max = small_number
		v_max = small_number
		w_max = small_number
		umax_loc = small_number
		vmax_loc = small_number
		wmax_loc = small_number
		umax_tmp = large_number

		ufmean_temploc = zero 
		umean_temploc = zero
		usmean_temploc = zero 

		mesh_veltemp = mesh_vel
		if (.not.movingcv) mesh_vel = zero
		if (I_AM_NODE_ZERO) write (*,'(A,3(2x,g17.8))') 'MESH VEL in compute_timestep: ', mesh_vel(:) 

		call calc_velreal(u, umean, ubcp)

		do k = 1, local_ni(3)
			do j = 1, local_ni(2)
				do i=1, local_ni(1)
					if (fluid_atijk(i,j,k)) then
						ufmean_temploc(:)= ufmean_temploc(:)+ ubcp(i,j,k,:)
					else
						usmean_temploc(:)= usmean_temploc(:) + ubcp(i,j,k,:)
					endif
					umean_temploc(:) = umean_temploc(:) + ubcp(i,j,k,:)

					umax_loc = MAX(umax_loc, ABS(ubcp(i,j,k,1)-mesh_vel(1)))
					vmax_loc = MAX(vmax_loc, ABS(ubcp(i,j,k,2)-mesh_vel(2)))
					wmax_loc = MAX(wmax_loc, ABS(ubcp(i,j,k,3)-mesh_vel(3)))
				enddo
			enddo
		enddo

		umax_loc = MAX(umax_loc, maxval(ABS(velbdy(:,1))))
		vmax_loc = MAX(vmax_loc, maxval(ABS(velbdy(:,2))))
		wmax_loc = MAX(wmax_loc, maxval(ABS(velbdy(:,3))))

		GLOBAL_DOUBLE_MAX(umax_loc,u_max,1,comm_cart_2d)
		GLOBAL_DOUBLE_MAX(vmax_loc,v_max,1,comm_cart_2d)
		GLOBAL_DOUBLE_MAX(wmax_loc,w_max,1,comm_cart_2d)

		umax_tmp = MAX(u_max,v_max)
		umax_tmp = MAX(umax_tmp,w_max)
		udiff = vis/dia_phys
		umax_tmp = MAX(umax_tmp,udiff)

		if (move_particles) then
			max_part_vel = zero

			partstart = 1
			do iphs = 1, nphases
				phase_mass(iphs) = rhos*pi*(phase_array(iphs)%dia)**3.d0/6.d0
				partend = partstart + phase_array(iphs)%npart- 1
				do m = partstart, partend
					if (is_mobile(m)) then
						part_vel_tmp(:) = force(m,:) / phase_mass(iphs) * dt + velbdy(m,:)

						do idim=1, ndim
							if ( abs(part_vel_tmp(idim)) > max_part_vel ) max_part_vel = abs(part_vel_tmp(idim))
						enddo
					endif
				enddo
				partstart = partend + 1
			enddo

			umax_tmp = MAX(umax_tmp, max_part_vel)
		endif

		!if (umax_tmp*dt/dx.gt.cfl) then
			dtnew = cfl*dx/umax_tmp
		!else
		!	dtnew = dt 
		!endif

		GLOBAL_DOUBLE_SUM(ufmean_temploc,ufmean_temp,3,comm_cart_2d)
		GLOBAL_DOUBLE_SUM(usmean_temploc,usmean_temp,3,comm_cart_2d)
		GLOBAL_DOUBLE_SUM(umean_temploc,umean_temp,3,comm_cart_2d)
		do idim = 1, ndim
			ufmean_temp(idim) = ufmean_temp(idim) / real(count_fluid,prcn)
			usmean_temp(idim) = usmean_temp(idim) / real(count_solid,prcn)
			umean_temp(idim)  = umean_temp(idim)  / real(global_n(1)*global_n(2)*global_n(3),prcn)
		enddo

		ufmean_old(:) = ufmean(:)
		ufmean(:) = ufmean_temp(:)
		dufmeandt = -cf*(ufmean(:)-ufmean_old(:))
		!dufmeandt = -cf*(ufmean_des(:)-ufmean(:))

		usmean = zero
		do iphs = 1, nphases
			phase_array(iphs)%mean_spec_vel(1:ndim) = zero
			partstart = phase_array(iphs)%pstart
			partend = phase_array(iphs)%pend
			do idim = 1, ndim
				phase_array(iphs)%mean_spec_vel(idim) = SUM(velbdy(partstart:partend,idim)) &
																	& /real(phase_array(iphs)%npart,prcn)
			enddo
			!!$ usmean(1:ndim) = usmean(1:ndim) + (phase_array(iphs)%volfracg)*(phase_array(iphs)%mean_spec_vel(1:ndim))
			usmean(1:ndim) = usmean(1:ndim) + (phase_array(iphs)%volfrac)*(phase_array(iphs)%mean_spec_vel(1:ndim))
			!write (*,*)'volfracg = ', (phase_array(iphs)%volfracg)
		enddo
		usmean(:) = usmean(:)/mean_volfrac
		mixmeanslip(:) = (one-maxvolfrac)*(usmean(:)-ufmean(:))
		mixmeanslipmod = DSQRT(doT_PRODUCT(mixmeanslip(1:ndim),mixmeanslip(1:ndim)))

		!do idim = 1, ndim
		!	if (nbody.gt.0) then
		!		usmean_act(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
		!	else
		!		usmean_act(idim) = zero
		!	endif
		!enddo
		usmean_act(:) = usmean(:)
    
		!usmean can become anything in IBM. SO we dont evaluate mean slip based on the current value of usmean.
		!As long as the particle velocities are correctly attained,
		!usmean_des has been reached and usmean calculated is irrelevant.!
		!new meanslip(:) = (ufmean(:)-usmean_des(:))*(one-maxvolfrac)
		meanslip(:) = (usmean_act(:)-ufmean(:))*(one-maxvolfrac)
		meanslipmod = DSQRT(dot_product(meanslip(1:ndim), meanslip(1:ndim)))
		re_out = mixmeanslipmod*char_length/vis

		if (I_AM_NODE_ZERO) then
			write (*,'(A25,3(2x,g17.8))')'USMEAN DES = ', usmean_des(:)
			write (*,'(A25,3(2x,g17.8))')'USMEAN ACT = ', usmean_act(:)
			write (*,'(A25,3(2x,g17.8))')'USMEAN MIX = ', usmean(:)
			write (*,'(A25,3(2x,g17.8))')'UFMEAN DES = ', ufmean_des(:)
			write (*,'(A25,3(2x,g17.8))')'UFMEAN ACTUAL = ', ufmean_temp(:)
			!write (*,'(A25,3(2x,g17.8))')'UFMEAN CC = ', ufmean_cc(:)
			write (*,'(A25,3(2x,g17.8))')'UMEAN ACTUAL = ', umean(:)
			write (*,'(A25,3(2x,g17.8))')'UMEAN TMP = ', umean_temp(:)
			write (*,'(A25,2(2x,g12.5))') "MAX VEL AND CFL:", umax_tmp,  umax_tmp*dt/dx
			write (*,'(A40,3(2x,g12.5))') "MEANSLIPMOD AND REY(MEANSLIPMOD):", mixmeanslipmod, re_out
		endif


		!Convert Pressure to real space
		call calc_pressure
    
		if (rks.eq.itrmax)then
			IF (adaptive.and..not.only_dem) then
				if (iglobstep>1) then
					dt_old = dt

					dt = dtnew*ab_relax + dt_old * (one-ab_relax)
					dt = min(dt,dtnew)
					!DT = DTNEW

					dt_ratio = dt/dt_old

					!FOR EULER TIME STEPPING 
					coef(1,1)=half
					coef(1,2)=half
					coef(1,3)=one+ dt_ratio/two
					coef(1,4)=-dt_ratio/two
					if (I_AM_NODE_ZERO) write (*,"(1a,1f)") "DT RATIO IN ADAPTIVE TIME STEPING = ", dt_ratio
				endif
			endif
		endif
		mesh_vel = mesh_veltemp
		if (iglobstep==1 .or. mod(iglobstep,skip_num)==0) then
			call reynolds_stress_tensor
			call compute_sijsij
		endif

		if (I_AM_NODE_ZERO) then
			call CPU_TIME (CPU1) 
			new_timestep_time_current = cpu1-cpu0
			new_timestep_time = new_timestep_time + new_timestep_time_current
		endif
	end subroutine compute_new_timestep



     
	subroutine interpolate_fields_to_new_mesh(delta_pos)
		use nlmainarrays, Only : ubcp, nlbcp, pbcp, onlbcp
		use field_tmp_arrays

		IMPLICIT NONE
		real(prcn), intent(in) :: delta_pos(ndim)
		integer :: i, j, k, vcell(ndim), pcell(ndim), ib,ie,jb,je,kb,ke, idim, onew
		real(prcn) :: new_pos(ndim), new_ppos(ndim), p_newmesh, u_newmesh(ndim), onll(ndim), nll(ndim), dfll(ndim), ppll(ndim),pmean

		do k = 1, mz
			new_pos(3) = k + delta_pos(3)
			new_ppos(3) = new_pos(3)
			do j = 1, my
				new_pos(2) = j + delta_pos(2)
				new_ppos(2) = new_pos(2)
				do i = 1, mx !new mx1
					new_pos(1) = i + delta_pos(1)
					new_ppos(1) = (i + 0.5 + delta_pos(1)) ! Position of the pressure grid point corresponding to vel grid point i
					new_ppos(1) = new_ppos(1) - 0.5

					do idim = 1, ndim
						if (new_pos(idim).lt.zero)then
							vcell(idim) = int(new_pos(idim)-1)
						else
							vcell(idim) = int(new_pos(idim))
						endif
						if (new_ppos(idim).lt.zero)then
							pcell(idim) = int(new_ppos(idim)-1)
						else
							pcell(idim) = int(new_ppos(idim))
						endif
					enddo
					p_newmesh = zero
					u_newmesh(1:ndim) = zero

					call interpolate_pdata(pcell,new_ppos,ppll,p_newmesh,1)
					call interpolate_udata(vcell,new_pos,ib,ie,jb,je,kb,ke,u_newmesh,nll,onll,dfll, 0,1, 1, onew) 
					nlbcp(i,j,k,1:ndim) = u_newmesh(1:ndim)
					onlbcp(i,j,k,1) = p_newmesh
!!$                if (iglobstep.eq.2)then
!!$                   write (*,*)'pbcp : ', pbcp(i,j,k)
!!$                   write (*,*)'p on newmesh :', p_newmesh
!!$                   stop
!!$                endif
				enddo
			enddo
		enddo


		!Store interpolated velocities and pressure fields into ubcp and pbcp
		do idim = 1, ndim
			do k = 1, mz
				do j = 1, my
					do i = 1, mx !new mx1
						ubcp(i,j,k,idim) = nlbcp(i,j,k,idim)
						if (idim.eq.1) pbcp(i,j,k) = onlbcp(i,j,k,idim)
					enddo
				enddo
			enddo
		enddo
		!new ubcp(mx,:,:,:) = ubcp(1,:,:,:)


#if 0
!new
		! Convert solutions back to complex space
		do idim = 1, ndim
			do i = 1, mx
				ur1(:,:) = ubcp(i,:,:,idim)-umean(idim)
				call ff2rc(ur1(:,:),u(i,:,:,idim))
				if (i.le.mx1) call ff2rc(pbcp(i,:,:),p(i,:,:))
			enddo
		enddo
#endif

		!fftw3 r2c
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					do i=1, mx
						urtmp(i,j,k) = ubcp(i,j,k,idim) - umean(idim)
					enddo
				enddo
			enddo
			call fftwr2c(urtmp, u(:,:,:,idim))
		enddo
		call fftwr2c(pbcp(1:local_ni(1),1:local_ni(2),1:local_ni(3)), p)
		!if (I_AM_NODE_ZERO) write (*,*) 'pmean = ', SUM(pbcp(:,:,:))/real(global_n(1)*global_n(2)*global_n(3),prcn)
	end subroutine interpolate_fields_to_new_mesh

  subroutine RUN_TIME_FILE_OPENER(funit, filename, formfile)
		use general_funcs 

		Implicit NOne

		Character(LEN=*) :: filename, formfile
		integer, intent(out) :: funit 

		logical :: filexist, isopen


		if (from_post) then
			OPEN(unit=funit,file="tmp_"//trim(filename), form=formfile,status='replace')
		else
			inQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)
			IF (.NOT.filexist) then
				funit = getnewunit(minunitno,maxunitno)
				OPEN(unit=funit,file=FILENAME, form=formfile,status='new')
			elseif (filexist.AND..NOT.isopen) then
				if (irestart.eq.0) then 
					funit = getnewunit(minunitno,maxunitno)
					OPEN(unit=funit,file=FILENAME, form=formfile,status='replace')
				elseif (irestart.eq.1) then 
					funit = getnewunit(minunitno,maxunitno)
					OPEN(unit=funit,file=FILENAME, form=formfile, POSITION="append")
				endif
			endif
		endif
  end subroutine RUN_TIME_FILE_OPENER

  subroutine RESTART_FILE_OPENER(funit, filename, formfile)
    use general_funcs 
    
    Implicit NOne 
    
    Character(LEN=*) :: filename, formfile
    integer, intent(out) :: funit 
    
    logical :: filexist, isopen
    
    inQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)
    
    IF (.NOT.filexist) then
       
       funit = getnewunit(minunitno,maxunitno)
       OPEN(unit=funit,file=FILENAME, form=formfile,status='new')
       
    elseif (filexist.AND..NOT.isopen) then
       funit = getnewunit(minunitno,maxunitno)
       
       OPEN(unit=funit,file=FILENAME, form=formfile,status='replace')
       
    endif
  end subroutine RESTART_FILE_OPENER


	subroutine calc_part_statistics(rks) 
		use general_funcs, only : instant_file_opener
		use init_turb, only : u_prime
		use general_funcs
		use randomno
		Implicit None 
		integer :: rks

		real(prcn) :: mean_drag(nphases), mean_force(nphases,ndim)
		real(prcn) :: fluctv(nbody,ndim),mean_vel(nphases,ndim),&
				& sdev_fluctv(nphases,ndim),sdev_fluctf(nphases,ndim),&
				& flupartslip(ndim), mixmeanvel(ndim), flupartslipmod,&
				& grant(nphases), phasicslip(nphsc2,ndim)
		real(prcn) :: mean_contact_force(ndim), fcvel_corr, char_force
		integer :: m,idim,iphs, partstart, partend, npart, phsc2count, jphs
		integer, SAVE :: velinfounit
		logical, SAVE :: FIRST_TIME = .true.
		CHARACTER*100 :: FILENAME, formfile

		if (.not.(rks.eq.itrmax))RETURN
#if 0
		if (I_AM_NODE_ZERO.and.first_time) then     
			FILENAME = TRIM(RUN_NAME)//'_part_info'//'.rst'
			formfile='unformatted'
			call  RUN_TIME_FILE_OPENER(unitpartinfo, FILENAME, formfile)

			velinfounit = getnewunit(minunitno,maxunitno)
			FILENAME = TRIM(RUN_NAME)//'_vel_info'//'.dat'
			formfile='formatted'
			call  RUN_TIME_FILE_OPENER(velinfounit, FILENAME, formfile)

			first_time = .false.
		endif
#endif

		filename = trim(run_name)//"_part_info"
		unitpartinfo = 1
		call instant_file_opener(filename,unitpartinfo,.false.)

		write (unitpartinfo) t-dt
		write (unitpartinfo) xc(1:nbody, 1:ndim)
		write (unitpartinfo) velbdy(1:nbody, 1:ndim)
		write (unitpartinfo) force(1:nbody, 1:ndim)
		write (unitpartinfo) pres(1:nbody, 1:ndim)
		write (unitpartinfo) visc(1:nbody, 1:ndim)
		write (unitpartinfo) contact_force(1:nbody, 1:ndim)
		write (unitpartinfo) frame_vel(1:ndim)
		write (unitpartinfo) frame_accln(1:ndim)
		write (unitpartinfo) ufmean(1:ndim)
		close (unitpartinfo)

		mixmeanvel = zero
		grant = zero
		do iphs = 1, nphases
			phase_array(iphs)%mean_spec_vel(1:ndim) = zero
			partstart = phase_array(iphs)%pstart
			partend = phase_array(iphs)%pend
			do idim = 1, ndim
				phase_array(iphs)%mean_spec_vel(idim) =&
				& SUM(velbdy(partstart:partend,idim))&
				&/real(phase_array(iphs)%npart,prcn)
			enddo
			mixmeanvel(1:ndim) = mixmeanvel(1:ndim) + (phase_array(iphs)&
				&%volfracg)*(phase_array(iphs)%mean_spec_vel(1:ndim))
		    
			do m = partstart, partend
				do idim = 1, ndim
					grant(iphs) = grant(iphs) + (phase_array(iphs)&
						&%mean_spec_vel(idim)-velbdy(m,idim))**2.d0
		   	enddo
			enddo
			grant(iphs) = grant(iphs)/(three*phase_array(iphs)%npart)
		enddo
		mixmeanvel(:) = mixmeanvel(:)/maxvolfrac
		flupartslip(:) = (usmean(:)-ufmean(:))
		flupartslipmod = doT_PRODUCT(flupartslip(1:ndim), flupartslip(1:ndim))
		flupartslipmod = DSQRT(flupartslipmod)
		 
		phsc2count = 0
		do iphs = 1, nphases-1
			do jphs = iphs+1, nphases
				phsc2count = phsc2count + 1
				phasicslip(phsc2count,1:ndim) = phase_array(jphs)&
						&%mean_spec_vel(1:ndim) - phase_array(iphs)&
						&%mean_spec_vel(1:ndim)
			enddo
		enddo

		do idim = 1, ndim
			mean_contact_force(idim) = SUM(contact_force(1:nbody,idim))/real(nbody,prcn)
		enddo

		if (zero_slip) then
			char_force = 3.d0*pi*vis*ucharmod*(one-maxvolfrac)*char_length
		else
			char_force = 3.d0*pi*vis*(mixmeanslipmod+small_number)*char_length
		endif

		fcvel_corr = zero
		do m = 1, nbody
			do idim = 1, ndim
				fcvel_corr = fcvel_corr + contact_force(m,idim)*(phase_array(iphs)&
				&%mean_spec_vel(idim)-velbdy(m,idim))
			enddo
		enddo
		fcvel_corr = fcvel_corr/real(nbody,prcn)


		filename = trim(run_name)//"_vel_info"
		velinfounit = 1
		call instant_file_opener(filename,velinfounit,.true.)

		if (zero_slip) then
			write (velinfounit,'(40(2x,g17.8))')t, t/t_conv, t/t_vis,&
				& (one-maxvolfrac)*ucharmod*char_length/vis,&
				& mixmeanvel(1:ndim)/ucharmod, ((phasicslip(iphs,idim)&
				&/ucharmod, idim = 1, ndim),iphs=1,phsc2count),&
				& DSQRT(grant(1:nphases))/ucharmod,&
				& mean_contact_force(1:ndim)/char_force, fcvel_corr&
				&/(char_force*ucharmod)
		else
			write (velinfounit,'(40(2x,g17.8))')t, t/t_conv, t/t_vis,&
				& (one-maxvolfrac)*flupartslipmod*char_length/vis,&
				& mixmeanvel(1:ndim)/flupartslipmod, ((phasicslip(iphs,idim)&
				&/flupartslipmod, idim = 1, ndim),iphs=1,phsc2count),&
				& DSQRT(grant(1:nphases))/flupartslipmod,&
				& mean_contact_force(1:ndim)/char_force, fcvel_corr&
				&/(char_force*flupartslipmod)
		endif
		close (velinfounit)
	end subroutine calc_part_statistics
  
subroutine gener_lattice_mod(nbody,domlin,xc,dia1, dbdy)

  use precision 
  use constants 
  use global_data, ONLY : RUN_NAME
  
  implicit none
  double precision, intent(in) ::domlin(3),dia1
  double precision, intent(out) :: xc(nbody,3)
  double precision ::dia, dmax, fac , ymin, ymax
  double precision :: rstep, rstep_buff, xtmp, ytmp, ztmp,lngth_spec1, dbdy(nbody), facdmax , dmaxbydmin
  integer, intent(in) :: nbody 
  integer :: nprob1 , i, j, k, ntot, ii, jj, kk, nx, ny, nz, np1, n
  real(prcn) :: ZP, XP, YP
  
  print*, 'in GENER LATTICE MOD, NBODY = ', NBODY
  
  !doml(:) = doMLin(:)! - MAXVAL(dbdy(1:nbody))
  dia = 1.2*dia1	! so that particles don't touch each other to begin with
  dmaxbydmin = MAXVAL(DBDY(1:nbody))/MinVAL(DBDY(1:nbody))
  if (ABS(dmaxbydmin-one).lt.small_number) then
     fac = 1.05
  else
     fac = 1.05
  endif
  

  !fac = 1.2
!!$    nx = (domlin(1)-half*dia)/dia
!!$    print*,'nx = ', nx,domlin(1)
!!$    nz = ceiling((real(doMLin(3)-half*dia)/dia)) 
!!$    print*,'nx = ', nz    
!!$    np1 = ceiling(real(nbody/nz))+ 1
!!$    ny = ceiling(real(np1/nx)) + 1
    
    !       Specifying the initial distribution of the particles
!!$    print*,'nx = ', nx, 'ny = ', ny, 'nz = ', nz
    n = 1
    dmax =  dia1!MAXVAL(dbdy(1:nbody))
    
    facdmax = 0.01*dmax

    yp = dmax*fac/two
    do While (n.lt.nbody) 
       !print*,'dmax = ', dmax 
       !nz = ceiling((real(doML(3)-dmax*fac)/(dmax*fac)))
       !nx = ceiling((doml(1)-dmax*fac)/(dmax*fac))
       nz = ceiling((real(doMLin(3)-dmax*fac)/(dmax*fac)))
       nx = ceiling((domlin(1)-dmax*fac)/(dmax*fac))
       

       do k = 1, nz
          zp =  dmax*half +  (k-1)*dmax*fac
          do i = 1, nx
             xp = dmax*half + (i-1)*dmax*fac
             XC(n,1) = xp
             XC(n,2) = yp
             XC(n,3) = zp
             if (N.Eq.NBODY) GOTO 200
             n = n+1
                       
          enddo
          
          !dmax =  MAXVAL(dbdy(n-1:nbody))
          !nx = floor((domlin(1)-dmax*fac)/(dmax*fac))
          !nx = ceiling((domlin(1)-half*dmax*fac)/dmax)

       enddo
       
       !dmax =  MAXVAL(dbdy(n-1:nbody))
       yp = yp + dmax*fac

200 CONTinUE       
    enddo
    



  open(2222, file = TRIM(RUN_NAME)//"_sphere_center_unfor.inp", form="unformatted", status="replace")

  open(unit=1000,file= TRIM(RUN_NAME)//"_lattice_dist.dat",form='formatted',status='replace')

  write (1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" '

  do i  = 1, nbody
     write (2222)xc(i,1),xc(i,2),xc(i,3),dbdy(i)*0.5d0
     
     write (1000,'(4(2x,f12.8))') xc(i,1),xc(i,2),xc(i,3),dbdy(i)/2.d0
  enddo

  close(2222, status = "keep")
  
  close(1000, status = "keep")
  
  
end subroutine gener_lattice_mod

  
subroutine gener_lattice(nbody1, nbody2,ibidisperse,volfrac1,volfrac2,dia1,dia2&
     &,percent_buff)

  use precision 
  use constants 
  use global_data, ONLY : RUN_NAME

  implicit none

  double precision, intent(in) :: percent_buff, dia1,  dia2, volfrac1, volfrac2
  double precision ::numdens1, numdens2
  double precision :: rstep, rstep_buff, xtmp, ytmp, ztmp,lngth_spec1
  integer, intent(out) :: nbody1, nbody2 
  integer :: nprob1 , i, j, k, ntot, ii, jj, kk

  logical, intent(in) ::ibidisperse
! dia2 = dia1*diaratio

  !    percent_buff = 1.0

  open(2222, file =  TRIM(RUN_NAME)//"_sphere_center_unfor.inp", form="unformatted", status="replace")

  open(unit=1000,file= TRIM(RUN_NAME)//"_lattice_dist.dat",form='formatted',status='replace')

  write (1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" '
!!!! bidisperse

  numdens1 = 6.*volfrac1*1./(pi*dia1**3)
  numdens2 = 6.*volfrac2*1./(pi*dia2**3)
  print*,'volffracs = ', volfrac1, volfrac2, dia1, dia2
  print*,' Numdens of first species..',numdens1
  print*,' Numdens of second species..',numdens2

!!!! since we assume the domain size to be unity (-0.5,+0.5)

  rstep  = dia1/2.

  !! the hard sphere code does not like
  !! touching spheres; so add a 1% of diameter buffer (0.5% for radius)

  rstep_buff = rstep + percent_buff/2.*rstep/100. 

  nprob1 = int(1./(2.*rstep_buff))
  write (*,*)'number of spheres along one direction...',nprob1

  write (*,*)'Number of rows needed..',numdens1/(nprob1*nprob1)	
  lngth_spec1 = (int(numdens1/(nprob1*nprob1))+1)*2.*rstep_buff

  write (*,*)'Total length along z required...',lngth_spec1

  xtmp = rstep_buff
  ytmp = rstep_buff
  ztmp = rstep_buff

  ntot = 0
  nbody1 = 0 
  nbody2 = 0
  do k = 1, nprob1
     do j = 1, nprob1
        do i = 1, nprob1

           write (2222)xtmp,ytmp,ztmp,dia1/2.d0

           write (1000,'(4(2x,f12.8))') xtmp-0.5,ytmp-0.5d0,ztmp-0.5d0,dia1/2.d0
           ntot = ntot + 1
           nbody1 = nbody1 + 1
           if (ntot.ge.numdens1) then 
              print*,'Total number of particles generated..',ntot
              ii = i
              jj = j
              kk = k
              print*,'final i, j, k...',ii,jj,kk
              goto 40
           endif

           xtmp  = xtmp + 2.*rstep_buff
        enddo
        xtmp = rstep_buff
        ytmp = ytmp + 2.*rstep_buff 
     enddo
     xtmp = rstep_buff
     ytmp = rstep_buff  
     ztmp = ztmp + 2.*rstep_buff 

     !if (ztmp.gt.volfrac1)then !this works because ztmp \in (0,1) and volfrac \in(0,1)

     if (ztmp.gt.lngth_spec1)then !this works because ztmp \in (0,1) and volfrac \in(0,1)
        print*,'!!!!!!! stop !!!!!!!!!! : first species has exceeded its domain'
        print*,' total number of particles generated...',ntot
        stop
     endif
  enddo

40 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! since we assume the domain size to be unity (-0.5,+0.5)
  if (ibidisperse)then

     rstep  = dia2/2.

     !! the hard sphere code does not like
     !! touching spheres; so add a 1% of diameter buffer (0.5% for radius)

     rstep_buff = rstep + percent_buff/2.*rstep/100. 

     nprob1 = int(1./(2.*rstep_buff))

     xtmp = rstep_buff
     ytmp = rstep_buff
     ztmp = rstep_buff

     ntot = 0

     do k = 1, nprob1
        do j = 1, nprob1
           do i = 1, nprob1
              write (2222)xtmp,ytmp,lngth_spec1+ztmp,dia2/2.d0

              nbody2 = nbody2 + 1

              write (1000,'(4(2x,f12.8))') xtmp-half,ytmp-half,lngth_spec1+ztmp-half,dia2/2.d0

              ntot = ntot + 1
              if (ntot.ge.numdens2) then 
                 print*,'Total number of particles generated..',ntot
                 goto 50
              endif

              xtmp  = xtmp + 2.*rstep_buff
           enddo
           xtmp = rstep_buff
           ytmp = ytmp + 2.*rstep_buff 
        enddo
        xtmp = rstep_buff
        ytmp = rstep_buff  
        ztmp = ztmp + 2.*rstep_buff 
        !if (volfrac1+ztmp.gt.1.0)then !this works because ztmp \in (0,1) and volfrac \in(0,1)
        if (lngth_spec1+ztmp.gt.1.0)then !this works because ztmp \in (0,1) and volfrac \in(0,1)
           print*,'!!!!!!! stop !!!!!!!!!! : second species has exceeded its domain'
           print*,' total number of particles generated...',ntot
           stop
	endif
     enddo

50   continue
  endif

  close(2222, status = "keep")

  close(1000, status = "keep")


!!!! volfrac of species based on number density

  print*,'Vol frac of first species...',pi*dia1**3/6.*nbody1,' input...',volfrac1
  print*,'Vol frac of second species...',pi*dia2**3/6.*nbody2, ' input...',volfrac2


end subroutine gener_lattice

subroutine scale_to_grid_units(nbody,npart,nphases,my,mbox,xperiodic,percent_buf&
     &,xp_in,rad_in, min_part_sep, toscale)

  use precision 
  use constants 
  use global_data, ONLY : RUN_NAME


  implicit none

  integer, intent(in) :: nphases
  integer, intent(inout) :: nbody
  integer, intent(inout),dimension(nphases):: npart
  integer, ALLOCATABLE, dimension(:) :: nsp
  logical, intent(in) ::  xperiodic
  logical, OPTIONAL :: toscale
  integer :: n,ndim ,i,idim,j,my,mbox,ixcount, partend,partstart, iphs

  real(prcn), intent(in),dimension(nphases) ::  percent_buf
  real(prcn), intent(in) :: min_part_sep
  real(prcn),intent(inout):: xp_in(nbody,3)
  real(prcn),intent(inout) :: rad_in(nbody)
  real(prcn) :: dist, overlap, tempx, tempy, tempz
  real(prcn) :: rmax(3), L(3), tempr(3), r
  real(prcn), dimension(nbody) :: radbdy
  real(prcn), dimension(nbody,3) :: xp
  integer, dimension(nbody)::  towrite	

  print*,'in SCALE TO GRID UNITS'
  !print*,'percent_buf = ', percent_buf(1:nphases)
  !nbody =SIZE(xp,1)
  ndim = SIZE(xp,2)
  xp(1:nbody,1:3) = xp_in(1:nbody,1:3)
  radbdy(1:nbody) = rad_in(1:nbody)

  rad_in(1:nbody) = zero 
  xp_in(1:nbody,1:3) = zero 

  ALLOCATE(nsp(nphases))

  n = nbody

  nsp(1:nphases) = npart(1:nphases)

  nbody = 0
  npart(1:nphases) = 0

  do i=1,n
     towrite (i)=1
  enddo

  partstart = 1
  if (.not.present(toscale))toscale = .true.
  
  if (toscale)then
     write (*,*)' SCALinG : ', toscale
     do iphs = 1, nphases
        partend = partstart + nsp(iphs)-1
        do i = partstart, partend
           do idim=1,ndim
              if (idim.eq.1)then
                 xp(i,idim)=xp(i,idim)*(mbox-1)+1.d0
              else
                 xp(i,idim)=xp(i,idim)*(my)+1
              endif
           enddo
           radbdy(i) = radbdy(i)/(one+percent_buf(iphs)*one/100.d0)
           radbdy(i) = radbdy(i)*(my)
        enddo
        partstart = partend + 1
     enddo
  endif

  write (*,*)'Max val of XP : ', MAXVAL(XP(1:n,1)), mbox, my
  write (*,*)'Max val of RADBDY : ', MAXVAL(RADBDY(1:n))
  !check if sphere centers are closer than min_part_sep grid points
  rmax(1) = (mbox-1)/two
  rmax(2:3) = my/two
  L(1) = mbox-1
  L(2:3) = my

  do i = 1, n
     do j = 1, n
        if (i.ne.j.and.(towrite (j).eq.1))then
           r = zero
           do idim = 1, ndim
              tempr(idim) = xp(i,idim) - xp(j,idim) ! compute the separation
              ! in each dimension
              if ((ABS(tempr(idim))).gt.rmax(idim)) then
                 if (tempr(idim).lt.zero) then
                    tempr(idim) = tempr(idim) + L(idim)
                 else
                    tempr(idim) = tempr(idim) - L(idim)
                 endif
              endif
              r = r + tempr(idim)**2.d0
           enddo
           r = DSQRT(r)
!!$           if (i.eq.2.and.j.eq.22)then
!!$              write (*,'(A20,3(2x,g12.5))')'Position of i = ', XP(I,1:3)
!!$              write (*,'(A20,3(2x,g12.5))')'Position of j = ', XP(J,1:3)
!!$              write (*,*) 'r = ', r, tempr(1:ndim), RMAX(1:NDIM)
!!$           endif
           if (r.le.radbdy(i)+radbdy(j)+min_part_sep)then
              towrite (i) = 0
              overlap = r - (radbdy(i)+radbdy(j)+min_part_sep)
              write (*,'(A20,2x,g12.5,A20, 2(2x,I6))')'Spheres closer than ',min_part_sep, '  grid points..', i,j
              write (*,'(A30,2x, g12.5)')'OVERLAP = ',  overlap/(radbdy(1))
              write (*,'(A20,3(2x,g12.5))')'Position of i = ', XP(I,1:3)
              write (*,'(A20,3(2x,g12.5))')'Position of j = ', XP(J,1:3)
           endif
        endif
     enddo
  enddo
  
  partstart = 1
  
  do iphs = 1, nphases
     partend = partstart + nsp(iphs)-1
     do i = partstart, partend
        if (towrite (i).eq.1)then
           nbody = nbody+1
           npart(iphs) = npart(iphs)+1
           rad_in(nbody) = radbdy(i)
           xp_in(nbody,1:3) = xp(i,1:3)
        endif
     enddo
     partstart = partend + 1
  enddo
  
end subroutine scale_to_grid_units

  end Module dependent_functions
