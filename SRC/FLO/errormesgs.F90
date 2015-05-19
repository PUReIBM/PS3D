!-------
! Dependent module
!-------
! Contains general error messages
!-------
! Code:   Steinli
! Author: Chidambaram Narayanan
!         Nuclear Engineering Laboratory
!         ETH Zurich
! Date:   June 2000
!-------
MODULE errormesgs
  Use precision         ! Independent module
  Use string_funcs
  Use general_funcs
  
  Use global_data       ! Dependent module

  Implicit none
  Private
  Public:: printerror
  
!--------
 CONTAINS
!--------

  !--------
  ! Centralized error messages
  !--------
  SUBROUTINE printerror(errortype,comment)
    Character(LEN=*):: errortype, comment

    SELECT CASE (errortype)
      CASE ("newunit")
        call separator(eunit,40,'e')
        write(eunit,*)'Error: ', comment
        write(eunit,*)'Unable to obtain unused unit number &
                     & to open file'
        call separator(eunit,40,'e')
      CASE DEFAULT
        call separator(eunit,40,'e')
        write(eunit,*)'Error: ', TRIM(errortype)
        write(eunit,*)TRIM(comment)
        call separator(eunit,40,'e')
    END SELECT
    STOP
  END SUBROUTINE printerror
END MODULE errormesgs

!-------
! List of routines
!-------
! o Subroutine printerror
!-------
