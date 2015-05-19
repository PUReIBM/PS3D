!---------
! Independent module
!---------
! o Contains all procedures pertaining to interpolation methods.
! o Uses Lagrangian Polynomial interpolation with options for
!   second, fourth, and sixth order interpolation.
! o Generic interface included one-dimensional and
!   three-dimensional possibilities. Two-d is not yet coded.
!-------
! Code:   Steinli
! Author: Chidambaram Narayanan
!         Nuclear Engineering Laboratory
!         ETH Zurich
! August 9, 2002
!-------
MODULE interpolation
  USE precision
  USE constants

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: interpolator, justweights

  !-------
  ! Generic interface
  !-------
  
  INTERFACE interpolator
    MODULE PROCEDURE interp_oned_scalar
    MODULE PROCEDURE interp_oned_vector
    MODULE PROCEDURE interp_threed_scalar
    MODULE PROCEDURE interp_threed_vector
    MODULE PROCEDURE interp_threed_scalar_yp
    MODULE PROCEDURE interp_threed_vector_yp
    MODULE PROCEDURE calc_weightderiv 
  END INTERFACE 

  !-------
  ! To avoid allocating for every call (saves time)
  ! Aug 8, 2002
  !-------
  SAVE
  INTEGER, PARAMETER:: maxorder=2
  REAL(prcn), DIMENSION(maxorder), TARGET:: xval, yval, zval
  REAL(prcn), DIMENSION(maxorder-1):: dx, dy, dz
  REAL(prcn), DIMENSION(maxorder,maxorder,maxorder)   &
            , TARGET:: weights

!-------
 CONTAINS
!-------
  
  !--------
  !--------
  ! Interpolate a scalar quantity in one space dimension.
  !--------
  !--------
  SUBROUTINE interp_oned_scalar(coor,scalar,ppos,interp_scl,weight_pointer)
    REAL(prcn), DIMENSION(:), INTENT(in):: coor, scalar
    REAL(prcn), INTENT(in):: ppos
    REAL(prcn), INTENT(out):: interp_scl
    REAL(prcn), DIMENSION(:), POINTER, OPTIONAL:: weight_pointer

    INTEGER:: order, iorig, i
    REAL(prcn):: zeta

    order = SIZE(coor)
    iorig = order/2

    DO i = 1,order-1
      dx(i) = coor(i+1) - coor(i)
    ENDDO
    zeta = (ppos - coor(iorig))/dx(iorig)

    !-------
    ! Get shape function values
    !-------
    SELECT CASE (order)
      CASE (2)
        DO i = 1,order
          xval(i) = shape2(zeta,i)
        ENDDO
      CASE (3)
        DO i = 1,order
          xval(i) = shape3(zeta,i,dx)
        ENDDO
      CASE (4)
        DO i = 1,order
          xval(i) = shape4(zeta,i,dx)
        ENDDO
      CASE (5)
        DO i = 1,order
          xval(i) = shape5(zeta,i,dx)
        ENDDO
      CASE (6)
        DO i = 1,order
          xval(i) = shape6(zeta,i,dx)
        ENDDO
    END SELECT

    !-------
    ! Calculate interpolated value
    !-------
    interp_scl = 0.0
    DO i = 1,order
      interp_scl = interp_scl + scalar(i)*xval(i)
    ENDDO

    !-------
    ! Return the weights (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
      weight_pointer => xval
    ENDIF
  END SUBROUTINE interp_oned_scalar

  !--------
  !--------
  ! Interpolate an arbitrary sized array in one space dimension.
  ! Uses the scalar interpolation to obtain the weights.
  !--------
  !--------
  SUBROUTINE interp_oned_vector(coor,vector,ppos,interp_vec,weight_pointer)
    REAL(prcn), DIMENSION(:), INTENT(in):: coor
    REAL(prcn), DIMENSION(:,:), INTENT(in):: vector
    REAL(prcn), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:), INTENT(out):: interp_vec
    REAL(prcn), DIMENSION(:), POINTER, OPTIONAL:: weight_pointer

    INTEGER:: order, vec_size, nv, i
    REAL(prcn), DIMENSION(:), POINTER:: weights_scalar

    order    = SIZE(coor)
    !print*,'In Interp_onedvector:order = ',order
    vec_size = SIZE(vector,2)

    !-------
    ! Interpolate first component and get weights
    !-------
    CALL interp_oned_scalar(coor,vector(:,1),ppos             &
                           ,interp_vec(1),weights_scalar)
 
    !-------
    ! Interpolate remaining components
    !-------
    DO nv = 2,vec_size
      interp_vec(nv) = 0.0
      DO i = 1,order
        interp_vec(nv) =  interp_vec(nv)  &
                        + vector(i,nv)*weights_scalar(i)
      ENDDO
    ENDDO

    !-------
    ! Return the weights (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
      weight_pointer => weights_scalar
    ENDIF
  END SUBROUTINE interp_oned_vector

  !--------
  !--------
  ! Interpolate a scalar quantity in three dimensions.
  !--------
  !--------
  SUBROUTINE interp_threed_scalar(coor,scalar,ppos,interp_scl,order&
       &,isch,weight_pointer) 
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor
    REAL(prcn), DIMENSION(:,:,:), INTENT(in):: scalar
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), INTENT(out):: interp_scl
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer
    CHARACTER*5, INTENT(in) :: isch
    REAL(prcn), DIMENSION(:,:), ALLOCATABLE:: zetacsi
    INTEGER:: i, j, k
    INTEGER, INTENT(in):: order
    INTEGER:: iorig
    REAL(prcn), DIMENSION(3):: zeta 

    !-------
    ! Get order of interpolation
    !-------
    ! 
    weights  = zero
    DO i = 1,order-1
       dx(i) = coor(i+1,1,1,1)-coor(i,1,1,1)
       dy(i) = coor(1,i+1,1,2)-coor(1,i,1,2)
       dz(i) = coor(1,1,i+1,3)-coor(1,1,i,3)
    ENDDO

    SELECT CASE(isch)

    CASE('lpi')
       !order = SIZE(coor,1)
       iorig = order/2

       !-------
       ! Find out center cell widths
       !-------

       !-------
       ! Zeta as defined in Bala/Maxey
       !-------
       zeta(1:3) = ppos(1:3) - coor(iorig,iorig,iorig,1:3)
       zeta(1) = zeta(1)/dx(iorig)
       zeta(2) = zeta(2)/dy(iorig)
       zeta(3) = zeta(3)/dz(iorig)

       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (2)
          DO i = 1,order
             xval(i) = shape2(zeta(1),i)
             yval(i) = shape2(zeta(2),i)
             zval(i) = shape2(zeta(3),i)
          ENDDO
       CASE (3)
          DO i = 1,order
             xval(i) = shape3(zeta(1),i,dx)
             yval(i) = shape3(zeta(2),i,dy)
             zval(i) = shape3(zeta(3),i,dz)
          ENDDO
       CASE (4)
          DO i = 1,order 
             ! print*, 'in interp....zetayp(3,i) = ', zetayp(3,i),zval(i),i
             xval(i) = shape4(zeta(1),i,dx)
             yval(i) = shape4(zeta(2),i,dy)
             zval(i) = shape4(zeta(3),i,dz)
             !Print*, ppos(1:3)
             !Print*,'int',i,xval(i),yval(i),zval(i)
!!$          xval(i) = shape4new(ppos(1),coor(1:order,1,1,1),i)
!!$          yval(i) = shape4new(ppos(2),coor(1,1:order,1,2),i)
!!$          zval(i) = shape4new(ppos(3),coor(1,1,1:order,3),i)
          ENDDO
       CASE (5)
          DO i = 1,order
             xval(i) = shape5(zeta(1),i,dx)
             yval(i) = shape5(zeta(2),i,dy)
             zval(i) = shape5(zeta(3),i,dz)
          ENDDO
       CASE (6)
          DO i = 1,order
             xval(i) = shape6(zeta(1),i,dx)
             yval(i) = shape6(zeta(2),i,dy)
             zval(i) = shape6(zeta(3),i,dz)
          ENDDO
       END SELECT

    CASE('csi')
      ! order = SIZE(coor,1)
       iorig = (order+1)/2

       !-------
       ! Find out center cell widths
       !-------
       ALLOCATE(zetacsi(3,order))
       !Zetacsi as defined in Yueng and Pope hence the name
       !The defintions for zetacsi are true only for a structured grid
       !if (order.eq.4) then 
      
       !end if
       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (4)
          DO i = 1, order
             zetacsi(1,i) = (-ppos(1) + coor(i,1,1,1))/dx(1)
             zetacsi(2,i) = (-ppos(2) + coor(1,i,1,2))/dy(1)
             zetacsi(3,i) = (-ppos(3) + coor(1,1,i,3))/dz(1)
          END DO
          DO i = 1,order
             
             xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
             yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
             zval(i) = shape4csi(zetacsi(3,i),i,dz,3) 
          ENDDO
       CASE(3)
          DO i = 1, order
             
             zetacsi(1,i) = ((-ppos(1) + coor(i,1,1,1))/dx(1))
             zetacsi(2,i) =((-ppos(2) + coor(1,i,1,2))/dy(1))
             zetacsi(3,i) = ((-ppos(3) + coor(1,1,i,3))/dz(1))
!!$             zetacsi(1,i) = (ppos(1) - coor(1,1,1,1))/(coor(order,1,1&
!!$                  &,1)-coor(1,1,1,1))
!!$             zetacsi(2,i) = (ppos(2) - coor(1,1,1,2))/(coor(1,order,1&
!!$                  &,2)-coor(1,1,1,2)) 
!!$             zetacsi(3,i) = (ppos(3) - coor(1,1,1,3))/(coor(1,1,order&
!!$                  &,3)-coor(1,1,1,3))
          END DO
          DO i = 1,order
             if((xval(1)-coor(1,1,1,1)).lt.dx(1)) then 
                xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
             else 
                xval(i) = shape3csiright(zetacsi(1,i),i,dx,1)
             endif
             if((yval(1)-coor(1,1,1,2)).lt.dy(1)) then 

                yval(i) = shape3csileft(zetacsi(2,i),i,dy,2) 
             else 
                
                yval(i) = shape3csiright(zetacsi(2,i),i,dy,2) 
             end if
             if((zval(1)-coor(1,1,1,3)).lt.dz(1)) then 
                zval(i) = shape3csileft(zetacsi(3,i),i,dz,3) 
             else
                
                zval(i) = shape3csiright(zetacsi(3,i),i,dz,3) 
             endif
             
             print*,'zeta = ',zetacsi(1,i), xval(i),i
          ENDDO
       END SELECT
       
       DEALLOCATE(zetacsi)
    END SELECT !SCHEME 
    !-------
    ! Calculate weights for the different nodes
    !-------
    DO 10 k = 1,order
    DO 10 j = 1,order
    DO 10 i = 1,order
      weights(i,j,k) = xval(i)*yval(j)*zval(k)
    10 CONTINUE
      !If(order.eq.3) Print*,'in interpo...sum wt=,',sum(weights),order

    !-------
    ! Calculate the interpolated value
    !-------
    interp_scl = 0.0
    DO 20 k = 1,order
    DO 20 j = 1,order
    DO 20 i = 1,order
      interp_scl = interp_scl + scalar(i,j,k)*weights(i,j,k)
    20 CONTINUE

    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
      weight_pointer => weights
    ENDIF
  END SUBROUTINE interp_threed_scalar

  !-------
  !-------
  ! Interpolate an arbitrary sized array in three dimensions.
  ! Uses the scalar interpolation to obtain the weights.
  !-------
  !-------
  SUBROUTINE interp_threed_vector(coor,vector,ppos,interp_vec,order,fun,weight_pointer) 
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor, vector
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:), INTENT(out):: interp_vec
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer
    CHARACTER*5 :: fun
    INTEGER, INTENT(in):: order
    INTEGER :: vec_size
    INTEGER:: i, j, k, nv
    REAL(prcn), DIMENSION(:,:,:), POINTER:: weights_scalar

    !-------
    ! Get order of interpolation
    !-------
    !order    = SIZE(coor,1)
    vec_size = SIZE(vector,4)
    !print*,'In Interp_threedvector:order = ',order,vec_size

    CALL interp_threed_scalar(coor,vector(:,:,:,1),ppos  &
                             ,interp_vec(1),order,fun,weights_scalar)

    !-------
    ! Calculate the interpolated velocities
    !-------
    DO nv = 2,vec_size
       interp_vec(nv) = 0.0
       DO  k = 1,order
          DO j = 1,order
             DO  i = 1,order
                interp_vec(nv) =  interp_vec(nv)  &
                     + vector(i,j,k,nv)*weights_scalar(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights_scalar
    ENDIF

  END SUBROUTINE interp_threed_vector
  
   SUBROUTINE interp_threed_scalar_yp(coor,scalar,ppos,interp_scl,weight_pointer,cs)
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor
    REAL(prcn), DIMENSION(:,:,:), INTENT(in):: scalar
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), INTENT(out):: interp_scl
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer
    CHARACTER(len=2) :: cs

    INTEGER:: i, j, k
    INTEGER:: order, iorig 
    REAL(prcn), DIMENSION(3,4):: zetayp !right now only for cubic spline interp rg 03/14/05

    !-------
    ! Get order of interpolation
    !-------
    order = SIZE(coor,1)
    iorig = order/2

    !-------
    ! Find out center cell widths
    !-------
    DO i = 1,order-1
      dx(i) = coor(i+1,1,1,1)-coor(i,1,1,1)
      dy(i) = coor(1,i+1,1,2)-coor(1,i,1,2)
      dz(i) = coor(1,1,i+1,3)-coor(1,1,i,3)
    ENDDO
    
    !Zetayp as defined in Yueng and Pope hence the name
    !The defintions for zetayp are true only for a structured grid
    IF (order.EQ.4) THEN 
       DO i = 1, order
          zetayp(1,i) = (-ppos(1) + coor(i,1,1,1))/dx(1)
          zetayp(2,i) = (-ppos(2) + coor(1,i,1,2))/dy(1)
          zetayp(3,i) = (-ppos(3) + coor(1,1,i,3))/dz(1)
       END DO
    END IF
    !-------
    ! Get shape function values
    !-------
    SELECT CASE (order)
    CASE (4)
       DO i = 1,order
          xval(i) = shape4csi(zetayp(1,i),i,dx,1)
          yval(i) = shape4csi(zetayp(2,i),i,dy,2)
          zval(i) = shape4csi(zetayp(3,i),i,dz,3) 
       ENDDO
    END SELECT

    !-------
    ! Calculate weights for the different nodes
    !-------
    DO 10 k = 1,order
    DO 10 j = 1,order
    DO 10 i = 1,order
      weights(i,j,k) = xval(i)*yval(j)*zval(k)
    10 CONTINUE

    !-------
    ! Calculate the interpolated value
    !-------
    interp_scl = 0.0
    DO 20 k = 1,order
    DO 20 j = 1,order
    DO 20 i = 1,order
      interp_scl = interp_scl + scalar(i,j,k)*weights(i,j,k)
    20 CONTINUE

    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
      weight_pointer => weights
    ENDIF
  END SUBROUTINE interp_threed_scalar_yp
  
  SUBROUTINE interp_threed_vector_yp(coor,vector,ppos,interp_vec,weight_pointer,cs)
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor, vector
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:), INTENT(out):: interp_vec
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer
    CHARACTER(len=2) :: cs
    INTEGER:: order, vec_size
    INTEGER:: i, j, k, nv
    REAL(prcn), DIMENSION(:,:,:), POINTER:: weights_scalar

    !-------
    ! Get order of interpolation
    !-------
    order    = SIZE(coor,1)
    vec_size = SIZE(vector,4)
    !print*,'In Interp_threedvector:order = ',order,vec_size

    CALL interp_threed_scalar_yp(coor,vector(:,:,:,1),ppos  &
                             ,interp_vec(1),weights_scalar,'cs')

    !-------
    ! Calculate the interpolated velocities
    !-------
    DO nv = 2,vec_size
       interp_vec(nv) = 0.0
       DO  k = 1,order
          DO j = 1,order
             DO  i = 1,order
                interp_vec(nv) =  interp_vec(nv)  &
                     + vector(i,j,k,nv)*weights_scalar(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights_scalar
    ENDIF

  END SUBROUTINE interp_threed_vector_yp
 
  SUBROUTINE calc_weightderiv(coor,ppos,weight_pointer,order, isch)
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor
    !Real(prcn), Dimension(:,:,:), Intent(in):: scalar
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    !Real(prcn), Intent(out):: interp_scl
    REAL(prcn), DIMENSION(:,:,:,:) :: weight_pointer
    CHARACTER(len=5), INTENT(in) :: isch
    INTEGER:: i, j, k,kk
    Real(prcn) :: dx1,dy1,dz1
    INTEGER, INTENT(in) :: order
    INTEGER :: iorig
    REAL(prcn), DIMENSION(3):: zeta
    REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: zetacsi !right now only
    ! for cubic spline interp rg 03/14/05 

    weight_pointer = zero
    !-------
    ! Get order of interpolation
    !-------
    !order = SIZE(coor,1)
    iorig = order/2
    ! print*,'in interp...Iorig = ',iorig !Debug

    !-------
    ! Find out center cell widths
    !-------
    DO i = 1,order-1
       dx(i) = coor(i+1,1,1,1)-coor(i,1,1,1)
       dy(i) = coor(1,i+1,1,2)-coor(1,i,1,2)
       dz(i) = coor(1,1,i+1,3)-coor(1,1,i,3)
    ENDDO

    !Zetacsi as defined in Yueng and Pope hence the name
    !The defintions for zetacsi are true only for a structured grid
    !if (order.eq.4) then 

    !end if
    !-------
    ! Zeta as defined in Bala/Maxey
    !-------


    !-------
    ! Get shape function values
    !-------
    SELECT CASE (isch)
    CASE ('csi') 
       Allocate(zetacsi(3,order))
      
       DO k = 1,3
          SELECT CASE(order)
          CASE(4)
             DO i = 1, order
                zetacsi(1,i) = (-ppos(1) + coor(i,1,1,1))/dx(1)
                zetacsi(2,i) = (-ppos(2) + coor(1,i,1,2))/dy(1)
                zetacsi(3,i) = (-ppos(3) + coor(1,1,i,3))/dz(1)
             END DO
             DO i = 1,order
                IF(k.EQ.1) THEN 
                   xval(i) = (shape4deriv_csi(zetacsi(1,i),i,dx))/dx(1)  
                   yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
                   zval(i) = shape4csi(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
                   yval(i) = (shape4deriv_csi(zetacsi(2,i),i,dy))/(dy(1))
                   zval(i) = shape4csi(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.3) THEN
                   xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
                   yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
                   zval(i) = (shape4deriv_csi(zetacsi(3,i),i,dz))/(dz(1))
                ENDIF
             ENDDO
          CASE(3)
             dx1 = coor(order,1,1,1)-coor(1,1,1,1)
             dy1 = coor(1,order,1,2)-coor(1,1,1,2)
             dz1 = coor(1,1,order,3)-coor(1,1,1,3)
             zetacsi(1,i) = (ppos(1) - coor(1,1,1,1))/dx1
             zetacsi(2,i) = (ppos(2) - coor(1,1,1,2))/dy1
             zetacsi(3,i) = (ppos(3) - coor(1,1,1,3))/dz1
             DO i = 1,order
                IF(k.EQ.1) THEN 
                   xval(i) = (shape3deriv_csi(zetacsi(1,i),i,dx))/dx1
                   yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
                   zval(i) = shape3csileft(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
                   yval(i) = (shape3deriv_csi(zetacsi(2,i),i,dy))/(dy1)
                   zval(i) = shape3csileft(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.3) THEN
                   xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
                   yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
                   zval(i) = (shape3deriv_csi(zetacsi(3,i),i,dz))/(dz1)
                ENDIF
             ENDDO
          END SELECT!order
          DO kk = 1,order
             DO j = 1,order
                DO i = 1,order
                   weight_pointer(i,j,kk,k) = xval(i)*yval(j)*zval(kk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO!end loop over the coordinate directions
       deallocate(zetacsi)
    CASE('lpi')
       zeta(1:3) = ppos(1:3) - coor(iorig,iorig,iorig,1:3)
       zeta(1) = zeta(1)/dx(iorig)
       zeta(2) = zeta(2)/dy(iorig)
       zeta(3) = zeta(3)/dz(iorig)
       DO k = 1,3
          DO i = 1,order
             IF(k.EQ.1) THEN 
                xval(i) = (shape4deriv(zeta(1),i,dx))/dx(iorig)  
                yval(i) = shape4(zeta(2),i,dy)
                zval(i) = shape4(zeta(3),i,dz)
             ELSEIF(k.EQ.2) THEN
                xval(i) = shape4(zeta(1),i,dx)
                yval(i) = (shape4deriv(zeta(2),i,dy))/(dy(iorig))
                zval(i) = shape4(zeta(3),i,dz)
             ELSEIF(k.EQ.3) THEN
                xval(i) = shape4(zeta(1),i,dx)
                yval(i) = shape4(zeta(2),i,dy)
                zval(i) = (shape4deriv(zeta(3),i,dz))/(dz(iorig))
             ENDIF
          ENDDO
          DO kk = 1,order
             DO j = 1,order
                DO i = 1,order
                   weight_pointer(i,j,kk,k) = -xval(i)*yval(j)*zval(kk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO!end loop over the coordinate directions
    END SELECT

  END SUBROUTINE calc_weightderiv

  !-------
  !-------
  ! To calculate just the weights using trilinear interpolation 
  !-------
  !-------
  FUNCTION justweights(coor,ppos)
    REAL(prcn), DIMENSION(:,:,:), POINTER:: justweights 
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(IN):: coor
    REAL(prcn), DIMENSION(3), INTENT(IN):: ppos

    INTEGER, PARAMETER:: order=2
    INTEGER:: i, j, k, iorig
    REAL(prcn):: dxl, dyl, dzl
    REAL(prcn), DIMENSION(3):: zeta

    iorig = order/2

    dxl = coor(order,1,1,1)-coor(order-1,1,1,1)
    dyl = coor(1,order,1,2)-coor(1,order-1,1,2)
    dzl = coor(1,1,order,3)-coor(1,1,order-1,3)

    zeta(1:3) = ppos(1:3) - coor(iorig,iorig,iorig,1:3)
    zeta(1) = zeta(1)/dxl
    zeta(2) = zeta(2)/dyl
    zeta(3) = zeta(3)/dzl

    DO i = 1,order
      xval(i) = shape2(zeta(1),i)
      yval(i) = shape2(zeta(2),i)
      zval(i) = shape2(zeta(3),i)
    ENDDO

    !-------
    ! Calculate weights for the different nodes
    !-------
    DO 10 k = 1,order
    DO 10 j = 1,order
    DO 10 i = 1,order
      weights(i,j,k) = xval(i)*yval(j)*zval(k)
    10 CONTINUE

    justweights => weights
  END FUNCTION justweights

  !-------
  !-------
  ! Second-order (linear) shape functions
  !-------
  !-------
  FUNCTION shape2(zeta,i)
    REAL(prcn):: shape2
    REAL(prcn):: zeta
    INTEGER:: i

    SELECT CASE (i)
      CASE (1)
        shape2 = one - zeta
      CASE (2)
        shape2 = zeta
    END SELECT
  END FUNCTION shape2

  !-------
  !-------
  ! Third-order (quadratic) shape functions
  !-------
  !-------
  FUNCTION shape3(zeta,i,dx)
    REAL(prcn):: shape3
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i

    REAL(prcn):: zh, num, denom

    SELECT CASE (i)
      CASE (1)
        zh    = dx(1)*zeta
        denom = dx(1)*(dx(1)+dx(2))
        num   = (zh - dx(1))*(zh - dx(1) - dx(2))
        shape3 = num/denom
      CASE (2)
        zh    = dx(1)*zeta
        denom = -dx(1)*dx(2)
        num   = zh*(zh - dx(1) - dx(2))
        shape3 = num/denom
      CASE (3)
        zh    = dx(1)*zeta
        denom = dx(2)*(dx(1)+dx(2))
        num   = zh*(zh - dx(1))
        shape3 = num/denom
    END SELECT
  END FUNCTION shape3

  !-------
  !-------
  ! Fourth-order (cubic) shape functions
  !-------
  !-------
  FUNCTION shape4(zeta,i,dx)
    REAL(prcn):: shape4
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i

    REAL(prcn):: r1, r2, c1, c2, c3

    r1 = dx(2)/dx(1)
    r2 = dx(3)/dx(1)

    SELECT CASE (i)
      CASE (1)
        c1 = 1.0/(1.0 + r1)
        c3 = 1.0/(1.0 + r1 + r1*r2)
        shape4 = -c1*c3*(r1**3.0)*(zeta)*(zeta-1.0)*(zeta-(1.0+r2))
        !shape4 = -(one/6.)*(zeta)*(zeta-one)*(zeta-two)
      CASE (2)
        c2 = 1.0/(1.0 + r2)
        shape4 = c2*(zeta-1.0)*(zeta*r1+1.0)*(zeta-(1.0+r2))
       ! shape4 = half*(zeta-one)*(zeta+one)*(zeta-two)
      CASE (3)
        c1 = 1.0/(1.0 + r1)
        shape4 = -(c1/r2)*(zeta)*(zeta*r1+1.0)*(zeta-(1.0+r2))
       ! shape4 = -half*(zeta)*(zeta+1)*(zeta-two)
      CASE (4)
        c2 = 1.0/(1.0 + r2)
        c3 = 1.0/(1.0 + r1 + r1*r2)
        shape4 = (c3*c2/r2)*(zeta)*(zeta*r1+1.0)*(zeta-1.0)
       ! shape4 = (one/6.)*(zeta)*(zeta+one)*(zeta-one)
    END SELECT
  END FUNCTION shape4

  !-------
  !-------
  ! Fifth-order (power of 4) shape functions
  !-------
  !-------
  FUNCTION shape5(zeta,i,dx)
    REAL(prcn):: shape5
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
   
    REAL(prcn):: d1, d2, d3, d4
    REAL(prcn):: num, denom, zh

    SELECT CASE (i)
      CASE (1)
        d1 = -dx(1)
        d2 = d1 - dx(2)
        d3 = d2 - dx(3)
        d4 = d3 - dx(4)
        denom = d1*d2*d3*d4
        zh = zeta*dx(2)
        num = zh*(zh -dx(2))*(zh -dx(2) -dx(3)) &
             *(zh -dx(2) -dx(3) -dx(4))
        shape5 = num/denom 
      CASE (2)
        d1 =  dx(1)
        d2 = -dx(2)
        d3 =  d2 - dx(3)
        d4 =  d3 - dx(4)
        denom = d1*d2*d3*d4
        zh = zeta*dx(2)
        num = (zh +dx(1))*(zh -dx(2))*(zh -dx(2) -dx(3)) &
             *(zh -dx(2) -dx(3) -dx(4))
        shape5 = num/denom 
      CASE (3)
        d1 =  dx(1) + dx(2)
        d2 =  dx(2)
        d3 = -dx(3)
        d4 =  d3 - dx(4)
        denom = d1*d2*d3*d4
        zh = zeta*dx(2)
        num = (zh +dx(1))*(zh)*(zh -dx(2) -dx(3)) &
             *(zh -dx(2) -dx(3) -dx(4))
        shape5 = num/denom 
      CASE (4)
        d1 =  dx(1) + dx(2) + dx(3)
        d2 =  d1 - dx(1)
        d3 =  d2 - dx(2)
        d4 = -dx(4)
        denom = d1*d2*d3*d4
        zh = zeta*dx(2)
        num = (zh +dx(1))*(zh)*(zh -dx(2)) &
             *(zh -dx(2) -dx(3) -dx(4))
        shape5 = num/denom 
      CASE (5)
        d1 =  dx(1) + dx(2) + dx(3) + dx(4)
        d2 =  d1 - dx(1)
        d3 =  d2 - dx(2)
        d4 =  d3 - dx(3)
        denom = d1*d2*d3*d4
        zh = zeta*dx(2)
        num = (zh +dx(1))*(zh)*(zh -dx(2)) &
             *(zh -dx(2) -dx(3))
        shape5 = num/denom 
    END SELECT
  END FUNCTION shape5

  !-------
  !-------
  ! Sixth-order (power of 5) shape functions
  !-------
  !-------
  FUNCTION shape6(zeta,i,dx)
    REAL(prcn):: shape6
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
   
    REAL(prcn):: d1, d2, d3, d4, d5
    REAL(prcn):: num, denom, zh

    SELECT CASE (i)
      CASE (1)
        d1 = -dx(1)
        d2 = d1 - dx(2)
        d3 = d2 - dx(3)
        d4 = d3 - dx(4)
        d5 = d4 - dx(5)
        denom = d1*d2*d3*d4*d5
        zh = zeta*dx(3)
        num = (zh + dx(2))*(zh)*(zh - dx(3)) &
             *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
        shape6 = num/denom 
      CASE (2)
        d1 =  dx(1)
        d2 = -dx(2)
        d3 =  d2 - dx(3)
        d4 =  d3 - dx(4)
        d5 =  d4 - dx(5)
        denom = d1*d2*d3*d4*d5
        zh = zeta*dx(3)
        num = (zh +dx(1) +dx(2))*(zh)*(zh -dx(3)) &
             *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
        shape6 = num/denom 
      CASE (3)
        d1 =  dx(1) + dx(2)
        d2 =  dx(2)
        d3 = -dx(3)
        d4 =  d3 - dx(4)
        d5 =  d4 - dx(5)
        denom = d1*d2*d3*d4*d5
        zh = zeta*dx(3)
        num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh -dx(3)) &
             *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
        shape6 = num/denom 
      CASE (4)
        d1 =  dx(1) + dx(2) + dx(3)
        d2 =  d1 - dx(1)
        d3 =  d2 - dx(2)
        d4 = -dx(4)
        d5 =  d4 - dx(5)
        denom = d1*d2*d3*d4*d5
        zh = zeta*dx(3)
        num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh) &
             *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
        shape6 = num/denom 
      CASE (5)
        d1 =  dx(1) + dx(2) + dx(3) + dx(4)
        d2 =  d1 - dx(1)
        d3 =  d2 - dx(2)
        d4 =  d3 - dx(3)
        d5 = -dx(5)
        denom = d1*d2*d3*d4*d5
        zh = zeta*dx(3)
        num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh) &
             *(zh -dx(3))*(zh - dx(3) -dx(4) -dx(5))
        shape6 = num/denom 
      CASE (6)
        d1 =  dx(1) + dx(2) + dx(3) + dx(4) + dx(5)
        d2 =  d1 - dx(1)
        d3 =  d2 - dx(2)
        d4 =  d3 - dx(3)
        d5 =  d4 - dx(4)
        denom = d1*d2*d3*d4*d5
        zh = zeta*dx(3)
        num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh) &
             *(zh -dx(3))*(zh - dx(3) -dx(4))
        shape6 = num/denom 
    END SELECT
  END FUNCTION shape6 

  FUNCTION shape3deriv_csi(zeta,i,dx)
    REAL(prcn):: shape3deriv_csi
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn) :: tmp
    
!!$    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN 
!!$       shape3deriv_csi = (half)*(two+zeta)**2.0
!!$    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN 
!!$       shape3deriv_csi = (one/six)*(-9.0*zeta**2.0-12.0*zeta)
!!$    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN 
!!$       shape3deriv_csi = (one/six)*(9.0*zeta**2.0-12.0*zeta)
!!$    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN 
!!$       shape3deriv_csi = (-half)*(two-zeta)**2.0
!!$    ELSE
!!$       shape3deriv_csi = zero
!!$    ENDIF 
    SELECT CASE (i)
    CASE (1)
       shape3deriv_csi = -two*(1-zeta)
    CASE (2)
       shape3deriv_csi = -two*zeta+two*(1-zeta)
    CASE (3)
       shape3deriv_csi = two*zeta
     END SELECT
   END FUNCTION shape3deriv_csi

  FUNCTION shape4deriv_csi(zeta,i,dx)
    REAL(prcn):: shape4deriv_csi
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn) :: tmp
    
    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN 
       shape4deriv_csi = (half)*(two+zeta)**2.0
    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN 
       shape4deriv_csi = (one/six)*(-9.0*zeta**2.0-12.0*zeta)
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN 
       shape4deriv_csi = (one/six)*(9.0*zeta**2.0-12.0*zeta)
    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN 
       shape4deriv_csi = (-half)*(two-zeta)**2.0
    ELSE
       shape4deriv_csi = zero
    ENDIF 
  END FUNCTION shape4deriv_csi

  FUNCTION shape4deriv(zeta,i,dx)
    REAL(prcn):: shape4deriv
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn) :: tmp
     
    REAL(prcn):: r1, r2, c1, c2, c3

    r1 = dx(2)/dx(1)
    r2 = dx(3)/dx(1)

    SELECT CASE (i)
      CASE (1)
        c1 = 1.0/(1.0 + r1)
        c3 = 1.0/(1.0 + r1 + r1*r2)
        tmp = -c1*c3*(r1**3.0)
        !shape4deriv = -c1*c3*(r1**3.0)*(zeta)*(zeta-1.0)*(zeta-(1.0+r2))
        shape4deriv = (1.0)*(zeta-1.0)*(zeta-(1.0+r2))&
             + (zeta)*(1.0)*(zeta-(1.0+r2))&
             +(zeta)*(zeta-1.0)*(1.0)
        shape4deriv = tmp*shape4deriv
      CASE (2)
        c2 = 1.0/(1.0 + r2)

        shape4deriv = (1.0)*(zeta*r1+1.0)*(zeta-(1.0+r2)) &
             +(zeta-1.0)*(r1)*(zeta-(1.0+r2)) &
             +(zeta-1.0)*(zeta*r1+1.0)*(1.0)
        shape4deriv = c2*shape4deriv
      CASE (3)
        c1 = 1.0/(1.0 + r1)
        tmp = -c1/r2
        shape4deriv = (1.0)*(zeta*r1+1.0)*(zeta-(1.0+r2))&
             +(zeta)*(r1)*(zeta-(1.0+r2))&
             +(zeta)*(zeta*r1+1.0)*(1.0)
        shape4deriv = tmp*shape4deriv
      CASE (4)
        c2 = 1.0/(1.0 + r2)
        c3 = 1.0/(1.0 + r1 + r1*r2)
        tmp = (c3*c2/r2)
        shape4deriv = (1.0)*(zeta*r1+1.0)*(zeta-1.0)&
             +(zeta)*(r1)*(zeta-1.0)&
             +(zeta)*(zeta*r1+1.0)*(1.0)
        shape4deriv = tmp*shape4deriv
    END SELECT
  END FUNCTION shape4deriv
  
  FUNCTION shape4csi(zeta,i,dx,dim)
    REAL(prcn):: shape4csi
    REAL(prcn),INTENT(in):: zeta
    REAL(prcn), DIMENSION(:),INTENT(in):: dx
    INTEGER,INTENT(in):: i,dim
    
    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN 
       shape4csi = (one/six)*(two+zeta)**3.0
    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN 
       shape4csi = (one/six)*(-3.0*zeta**3.0-six*zeta**2.0+4.0)
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN 
       shape4csi = (one/six)*(3.0*zeta**3.0-six*zeta**2.0+4.0)
    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN 
       shape4csi = (one/six)*(two-zeta)**3.0
    ELSE
       shape4csi = zero
       !print*,'shape4rg .... zeta=',zeta,dim
    ENDIF
   
  END FUNCTION shape4csi 

  FUNCTION shape3csileft(zeta,i,dx,dim)
    REAL(prcn):: shape3csileft
    REAL(prcn),INTENT(in):: zeta
    REAL(prcn), DIMENSION(:),INTENT(in):: dx
    INTEGER,INTENT(in):: i,dim
    
    IF (zeta.GE.-one.AND.zeta.LE.zero) THEN 
       shape3csileft = (half)*(zeta+one)**2.0
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN 
       shape3csileft = (half)*(-two*zeta**2.+2.0*zeta+1.0)
    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN 
       shape3csileft = (half)*(zeta-two)**2.0
    ELSE
       shape3csileft = zero
       !print*,'shape4rg .... zeta=',zeta,dim
    ENDIF

!!$    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN 
!!$       shape3csileft = (half)*(two+zeta)**2.0
!!$    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN 
!!$       shape3csileft = (half)*(-two*zeta**2.-two*zeta+one)
!!$    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN 
!!$       shape3csileft = (half)*(one-zeta)**2.0
!!$    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN 
!!$       shape3csileft = (half)*(two-zeta)**2.0
!!$    ELSE
!!$       shape3csileft = zero
!!$       !print*,'shape4rg .... zeta=',zeta,dim
!!$    ENDIF
!!$    SELECT CASE (i)
!!$    CASE (1)
!!$       shape3csileft = (1-zeta)**2.
!!$       !shape4 = -(one/6.)*(zeta)*(zeta-one)*(zeta-two)
!!$    CASE (2)
!!$        shape3csileft = two*zeta*(1-zeta)
!!$     CASE (3)
!!$        shape3csileft = zeta**2.
!!$     END SELECT
  END FUNCTION shape3csileft


  FUNCTION shape3csiright(zeta,i,dx,dim)
    REAL(prcn):: shape3csiright
    REAL(prcn),INTENT(in):: zeta
    REAL(prcn), DIMENSION(:),INTENT(in):: dx
    INTEGER,INTENT(in):: i,dim
    
  

    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN 
       shape3csiright = (half)*(two+zeta)**2.0
    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN 
       shape3csiright = (half)*(-two*zeta**2.-two*zeta+one)
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN 
       shape3csiright = (half)*(one-zeta)**2.0
    ELSE
       shape3csiright = zero
       !print*,'shape4rg .... zeta=',zeta,dim
    ENDIF
  END FUNCTION shape3csiright

  FUNCTION shape4new(pos,x,i)
    REAL(prcn):: shape4new
    REAL(prcn):: pos
    REAL(prcn), DIMENSION(:):: x
    INTEGER:: i
    
    REAL(prcn):: r1, r2,num,den
     
    SELECT CASE (i)
    CASE (1) 
       num = (pos-x(2))*(pos - x(3))*(pos - x(4))
       den = (x(1) - x(2))*(x(1) - x(3))*(x(1) - x(4))
       shape4new = num/den
    CASE (2)
       num = (pos-x(1))*(pos - x(3))*(pos - x(4))
       den = (x(2) - x(1))*(x(2) - x(3))*(x(2) - x(4))
       shape4new = num/den 
    CASE (3)
       num = (pos-x(1))*(pos - x(2))*(pos - x(4))
       den = (x(3) - x(1))*(x(3) - x(2))*(x(3) - x(4))
       shape4new = num/den 
    CASE (4)
       num = (pos-x(1))*(pos - x(2))*(pos - x(3))
       den = (x(4) - x(1))*(x(4) - x(2))*(x(4) - x(3))
       shape4new = num/den 
    END SELECT
  END FUNCTION shape4new
END MODULE interpolation


!---------
! List of functions
!---------
! o subroutine interpolator
! o function justweights
! o function shape2
! o function shape4
! o function shape6
!---------
