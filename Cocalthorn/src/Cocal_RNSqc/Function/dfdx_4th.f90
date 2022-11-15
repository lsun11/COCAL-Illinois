function dfdx_4th(x,f,v)
  use phys_constant, only : long
  implicit none
  real(long) :: dfdx_4th
  real(long) :: x(5),f(5), v
  real(long) :: dx01, dx02, dx03, dx04
  real(long) :: dx10, dx12, dx13, dx14
  real(long) :: dx20, dx21, dx23, dx24
  real(long) :: dx30, dx31, dx32, dx34
  real(long) :: dx40, dx41, dx42, dx43
  real(long) ::  xv0,  xv1,  xv2,  xv3,  xv4
  real(long) :: wex0, wex1, wex2, wex3, wex4
  integer    ::  ir0,  ir1,  ir2,  ir3,  ir4
!
  ir0 = 1
  ir1 = 2
  ir2 = 3
  ir3 = 4
  ir4 = 5
  dx01 = x(ir0) - x(ir1)
  dx02 = x(ir0) - x(ir2)
  dx03 = x(ir0) - x(ir3)
  dx04 = x(ir0) - x(ir4)
  dx12 = x(ir1) - x(ir2)
  dx13 = x(ir1) - x(ir3)
  dx14 = x(ir1) - x(ir4)
  dx23 = x(ir2) - x(ir3)
  dx24 = x(ir2) - x(ir4)
  dx34 = x(ir3) - x(ir4)
  dx10 = - dx01
  dx20 = - dx02
  dx21 = - dx12
  dx30 = - dx03
  dx31 = - dx13
  dx32 = - dx23
  dx40 = - dx04
  dx41 = - dx14
  dx42 = - dx24
  dx43 = - dx34
  xv0 = v - x(ir0)
  xv1 = v - x(ir1)
  xv2 = v - x(ir2)
  xv3 = v - x(ir3)
  xv4 = v - x(ir4)
!
  wex0 = (xv2*xv3*xv4 + xv1*xv3*xv4 + &
  &       xv1*xv2*xv4 + xv1*xv2*xv3)/(dx01*dx02*dx03*dx04)
  wex1 = (xv2*xv3*xv4 + xv0*xv3*xv4 + &
  &       xv0*xv2*xv4 + xv0*xv2*xv3)/(dx10*dx12*dx13*dx14)
  wex2 = (xv1*xv3*xv4 + xv0*xv3*xv4 + &
  &       xv0*xv1*xv4 + xv0*xv1*xv3)/(dx20*dx21*dx23*dx24)
  wex3 = (xv1*xv2*xv4 + xv0*xv2*xv4 + &
  &       xv0*xv1*xv4 + xv0*xv1*xv2)/(dx30*dx31*dx32*dx34)
  wex4 = (xv1*xv2*xv3 + xv0*xv2*xv3 + &
  &       xv0*xv1*xv3 + xv0*xv1*xv2)/(dx40*dx41*dx42*dx43)
  dfdx_4th = wex0*f(ir0) + wex1*f(ir1) + wex2*f(ir2) & 
  &        + wex3*f(ir3) + wex4*f(ir4)
!
end function dfdx_4th
