Module nlphimainarrays
  Use precision 
  Use constants 
  Use global_data
  implicit none
  !Save 
  Public
   
  real(prcn), DImension(:,:,:,:),pointer ::  onlbcp,nlbcp, ubcp

!!$  real(prcn) ::  ubc(mxf,my,mz,ndim)
!!$  real(prcn) ::  nlbc(mxf,my,mz,ndim)
!!$  real(prcn) ::  onlbc(mxf,my,mz,ndim)
!!$  real(prcn) ::  pbc(mxf,my,mz)
end Module nlmainarrays
