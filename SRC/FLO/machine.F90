MODULE machine
  USE global_data 

  !              record length used in open statement for unformatted,
  !              direct access file, with 512 bytes per record
  INTEGER  OPEN_N1
  !
  !              number of DOUBLE PRECISION words in 512 bytes
  INTEGER  NWORDS_DP
  !
  !              number of REAL words in 512 bytes
  INTEGER  NWORDS_R
  !
  !              number of INTEGER words in 512 bytes
  INTEGER  NWORDS_I
  !
  LOGICAL :: JUST_FLUSH = .TRUE.

CONTAINS

  !
  SUBROUTINE MACHINE_CONS
    !
    !
    IMPLICIT NONE
    !
    OPEN_N1   = 512
    NWORDS_DP =  64
    NWORDS_R  = 128
    NWORDS_I  = 128
    JUST_FLUSH = .TRUE.
    !
    RETURN
  END SUBROUTINE MACHINE_CONS
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !                                                                      C
  !  Module name: GET_RUN_ID                                             C
  !  Purpose: get the run id for this run                                C
  !                                                                      C
  !  Author: P. Nicoletti                               Date: 16-DEC-91  C
  !  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
  !                                                                      C
  !  Revision Number: 1                                                  C
  !  Purpose: add ndoe name                                              C
  !  Author: P.Nicoletti                                Date: 07-FEB-92  C
  !  Reviewer:                                          Date: dd-mmm-yy  C
  !                                                                      C
  !  Literature/Document References:                                     C
  !                                                                      C
  !  Variables referenced: None                                          C
  !  Variables modified: ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, ID_MINUTE   C
  !                      ID_SECOND, ID_NODE                              C
  !                                                                      C
  !  Local variables: TIME_ARRAY                                         C
  !                                                                      C
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  !
  SUBROUTINE GET_RUN_ID
    !
    IMPLICIT NONE
    !
    !             temporary array to hold time data
    INTEGER DAT(8)
    CHARACTER*10 DATE, TIM, ZONE


    ! Intel Linux compiler supports this function thru it's portability library
    CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
    ID_YEAR   = DAT(1)
    ID_MONTH  = DAT(2)
    ID_DAY    = DAT(3)
    ID_HOUR   = DAT(5)
    ID_MINUTE = DAT(6)
    ID_SECOND = DAT(7)



    ! Intel Linux compiler supports this function thru it's portability library
    call hostnm(ID_NODE)      
    !
    RETURN
  END SUBROUTINE GET_RUN_ID
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !                                                                      C
  !  Module name: CPU_TIME (CPU)                                         C
  !  Purpose: get the CPU time for the run                               C
  !                                                                      C
  !  Variables referenced: None                                          C
  !  Variables modified: None                                            C
  !                                                                      C
  !  Local variables: TA, XT                                             C
  !                                                                      C
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  !
  SUBROUTINE CPU_TIME(CPU)
    !
    IMPLICIT NONE
    !
    ! passed arguments
    !
    !                      cpu time since start of run
    DOUBLE PRECISION CPU

    INTEGER, SAVE :: COUNT_OLD=0, WRAP=0
    !
    ! local variables
    !

    !                       clock cycle
    INTEGER           COUNT

    !                       number of cycles per second
    INTEGER           COUNT_RATE

    !                       max number of cycles, after which count is reset to 0
    INTEGER           COUNT_MAX

    !
    ! Intel Linux compiler supports this function thru it's portability library
    CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
    IF(COUNT_OLD .GT. COUNT) THEN
       WRAP = WRAP + 1
    ENDIF
    COUNT_OLD = COUNT

    CPU           = DBLE(COUNT)/DBLE(COUNT_RATE) &
         + DBLE(WRAP) * DBLE(COUNT_MAX)/DBLE(COUNT_RATE)


    RETURN
  END SUBROUTINE CPU_TIME
end MODULE machine
