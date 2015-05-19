MODULE meanflo
#include "ibm.h"
  USE precision 
  USE constants 
  USE global_data
  USE tridiagonal
  USE gauss_seidel
  IMPLICIT NONE 
CONTAINS

  SUBROUTINE calc_meanflo(rks)
    
    !       this is a 3d code
    
    IMPLICIT NONE 
    INTEGER, INTENT(in) ::  rks
    
    !-----------------------------------------------------------------------
    !	local variables
    
    INTEGER :: j
    
    umean(1) = umean(1) + (-mpg(1)+frmean(1))*dt
    umean(2) = umean(2) + (-mpg(2)+frmean(2))*dt
    umean(3) = umean(3) + (-mpg(3)+frmean(3))*dt
    
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A25,3(2x,g12.5))')'FR MEAN   = ', (FRMEAN(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g12.5))')'MPG   = ', (MPG(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g12.5))')'UMEAN = ', (UMEAN(j), j = 1, ndim)
    end if
    RETURN
  end SUBROUTINE calc_meanflo
end MODULE meanflo


    
