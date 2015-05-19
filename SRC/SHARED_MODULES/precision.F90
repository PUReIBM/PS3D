!-------
! Independent module
!-------
! To set the precision for the whole code
!-------
! Code:   Steinli
! Author: Chidambaram Narayanan
!         Nuclear Engineering Laboratory
!         ETH Zurich
!-------
MODULE precision
  Implicit none
  Save
  Public
  Integer, Parameter:: prcn=8
  Integer, Parameter:: iprec=8
END MODULE precision
