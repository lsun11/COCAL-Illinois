function dfdx_2nd(x,f,v)
  implicit none
  real(8) :: dfdx_2nd
  real(8) :: x(3),f(3), v
  real(8) :: dx01, dx02, dx12
  real(8) :: dx10, dx20, dx21
  real(8) ::  xv0,  xv1,  xv2
  real(8) :: wex0, wex1, wex2
  integer ::  ir0,  ir1,  ir2
!
  ir0 = 1
  ir1 = 2
  ir2 = 3
  dx01 = x(ir0) - x(ir1)
  dx02 = x(ir0) - x(ir2)
  dx12 = x(ir1) - x(ir2)
  dx10 = - dx01
  dx20 = - dx02
  dx21 = - dx12
  xv0 = v - x(ir0)
  xv1 = v - x(ir1)
  xv2 = v - x(ir2)
!
  wex0 = (xv2 + xv1)/(dx01*dx02)
  wex1 = (xv2 + xv0)/(dx10*dx12)
  wex2 = (xv1 + xv0)/(dx20*dx21)
  dfdx_2nd = wex0*f(ir0) + wex1*f(ir1) + wex2*f(ir2)
!
end function dfdx_2nd
