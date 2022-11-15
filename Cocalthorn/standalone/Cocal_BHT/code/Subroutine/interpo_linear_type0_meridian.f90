! linear interpolation on the GR coordinate.
subroutine interpo_linear_type0_meridian(val,fnc,ir,it,ip)
  use phys_constant, only : long
  implicit none
  real(long) :: val
  real(long), pointer :: fnc(:,:,:)
  integer :: ir, it, ip
  val = 0.25d0*(fnc(ir,it  ,ip) + fnc(ir-1,it  ,ip)  &
  &           + fnc(ir,it-1,ip) + fnc(ir-1,it-1,ip))
end subroutine interpo_linear_type0_meridian
