!-------
! Independent module
!-------
! Code:   Steinli
! Author: Chidambaram Narayanan
!         Nuclear Engineering Laboratory
!         ETH Zurich
!-------
MODULE constants
  USE precision

  Implicit none
  Save
  Public

  Real(iprec), Parameter:: zero = 0.0_iprec
  Real(iprec), Parameter:: one  = 1.0_iprec
  Real(iprec), Parameter:: half = 0.5_iprec
  Real(iprec), Parameter:: two = 2.0_iprec
  
  Real(iprec), Parameter:: three = 3.0_iprec
  Real(iprec), Parameter:: four = 4.0_iprec
  Real(iprec), Parameter:: five = 4.0_iprec
  Real(iprec), Parameter:: six = 6.0_iprec
  Real(iprec), Parameter:: seven = 7.0_iprec
  Real(iprec), Parameter:: eight = 8.0_iprec
  Real(iprec), Parameter:: nine = 9.0_iprec
  Real(iprec), Parameter:: ten = 10.0_iprec
  
  real(iprec), Parameter ::  qrt=0.25_iprec
  real(iprec), Parameter :: thrd = (1.0_iprec/3.0_iprec)
  real(iprec), Parameter :: thqrt = 0.75_iprec
  Real(iprec), Parameter:: gravity=9.81_iprec

  Real(iprec), Parameter:: tenm2=1.0e-2_iprec
  Real(iprec), Parameter:: tenm3=1.0e-3_iprec
  Real(iprec), Parameter:: tenm4=1.0e-4_iprec
  Real(iprec), Parameter:: tenm5=1.0e-5_iprec
  Real(iprec), Parameter:: tenm6=1.0e-6_iprec
  Real(iprec), Parameter:: tenm7=1.0e-7_iprec
  Real(iprec), Parameter:: tenm8=1.0e-8_iprec
  Real(iprec), Parameter:: tenm9=1.0e-9_iprec
  Real(iprec), Parameter:: infty=1.0e+35_iprec
  Real(iprec), Parameter:: negli=1.0e-35_iprec

  Real(iprec), Parameter:: milli=1.0e-3_iprec
  Real(iprec), Parameter:: micro=1.0e-6_iprec
  Real(iprec), Parameter:: nano=1.0e-9_iprec
  Real(iprec), Parameter:: angstrom=1.0e-10_iprec

  Real(iprec):: pi, twopi, onebypi, fourthirdpi

  Complex(iprec), Parameter:: ei=(zero,one)
  Complex(iprec), Parameter:: czero=(zero,zero)
  Complex(iprec), Parameter:: cone=(one,zero)

  REAL(PRCN), PARAMETER :: SMALL_NUMBER = 1.0E-35
  REAL(PRCN), PARAMETER :: LARGE_NUMBER = 1.0E+32

  REAL(PRCN), PARAMETER :: oneeighty = 180.d0

  INTEGER, PARAMETER :: UNDEFINED_I = 12345678
  REAL(PRCN), PARAMETER :: UNDEFINED_R = 12345678.D0

!-------
 CONTAINS
!-------

  SUBROUTINE calculate_constants
    Implicit none

    pi = 4.0_iprec*ATAN(one)
    twopi = 2.0*pi
    onebypi = 1.0_iprec/pi
    fourthirdpi = (four/three)*pi
  END SUBROUTINE calculate_constants
END MODULE constants
