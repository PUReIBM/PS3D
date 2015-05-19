module reversal_pts
  USE PRECISION
  
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC:: revpt, parrayrev, revpt_position, update_cell_index
  
  INTEGER, PARAMETER, Save:: maxorder=6
  
  !object reversal point
  Type revpt
     Private
     REAL(prcn), DIMENSION(3):: position
     INTEGER, DIMENSION(3):: cell
     REAL(prcn), DIMENSION(maxorder,maxorder,maxorder) :: wts
  END TYPE revpt
  
  !-------
  ! Have to define a new type (reversal point-array pointer)
  ! in order to declare an array of reversal points-array pointers
  !-------
  TYPE parrayrev
     TYPE(revpt), DIMENSION(:), POINTER:: p
  END TYPE parrayrev
  
  
  !-------
  ! Generic interface definitions
  !-------
  INTERFACE revpt_position
     MODULE PROCEDURE revpt_position_scalar
     MODULE PROCEDURE revpt_position_vector
  END INTERFACE
  
  !-------
  !-------
  ! Access the reversal point position vector
  !-------
  !-------
Contains 
  FUNCTION revpt_position_vector(p)
    REAL(prcn), DIMENSION(3):: revpt_position_vector
    TYPE(revpt), INTENT(IN)::p
    
    revpt_position_vector = p%position
  END FUNCTION revpt_position_vector
  
  !-------
  !-------
  ! Access the individual reversal point  position
  !-------
  !-------
  FUNCTION revpt_position_scalar(p,idirection)
    REAL(prcn):: revpt_position_scalar
    TYPE(revpt), INTENT(IN)::p
    INTEGER:: idirection

    revpt_position_scalar = p%position(idirection)
  END FUNCTION revpt_position_scalar

  !subroutine to update the rev pt cell index. 
  SUBROUTINE update_cell_index(p,newcell)
    TYPE(revpt), INTENT(INOUT):: p
    INTEGER, DIMENSION(3):: newcell
    
    p%cell = newcell
  END SUBROUTINE update_cell_index
end module reversal_pts
