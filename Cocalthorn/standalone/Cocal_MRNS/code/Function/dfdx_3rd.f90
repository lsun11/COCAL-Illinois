function dfdx_3rd(x,f,v)
  implicit none
  real(8) :: dfdx_3rd
  real(8) :: x(4),f(4), v
  real(8) :: dx01, dx02, dx03, dx10, dx12, dx13
  real(8) :: dx20, dx21, dx23, dx30, dx31, dx32
  real(8) :: xv0, xv1, xv2, xv3
  real(8) :: wex0, wex1, wex2, wex3
  integer :: ir0, ir1, ir2, ir3
!
  ir0 = 1
  ir1 = 2
  ir2 = 3
  ir3 = 4
  dx01 = x(ir0) - x(ir1)
  dx02 = x(ir0) - x(ir2)
  dx03 = x(ir0) - x(ir3)
  dx12 = x(ir1) - x(ir2)
  dx13 = x(ir1) - x(ir3)
  dx23 = x(ir2) - x(ir3)
  dx10 = - dx01
  dx20 = - dx02
  dx21 = - dx12
  dx30 = - dx03
  dx31 = - dx13
  dx32 = - dx23
  xv0 = v - x(ir0)
  xv1 = v - x(ir1)
  xv2 = v - x(ir2)
  xv3 = v - x(ir3)
  wex0 = (xv2*xv3 + xv1*xv3 + xv1*xv2)/(dx01*dx02*dx03)
  wex1 = (xv2*xv3 + xv0*xv3 + xv0*xv2)/(dx10*dx12*dx13) 
  wex2 = (xv1*xv3 + xv0*xv3 + xv0*xv1)/(dx20*dx21*dx23) 
  wex3 = (xv1*xv2 + xv0*xv2 + xv0*xv1)/(dx30*dx31*dx32) 
  dfdx_3rd = wex0*f(ir0) + wex1*f(ir1) + wex2*f(ir2) + wex3*f(ir3)
!
end function dfdx_3rd
