subroutine interpo_linear_surface_type0(val,fnc,ir,it,ip)
  use phys_constant, only : long
  implicit none
  real(long) :: val
  real(long), pointer :: fnc(:,:,:)
  integer :: ir, it, ip
  val = 0.25d0*(fnc(ir,it,ip)   + fnc(ir,it-1,ip)  &
              + fnc(ir,it,ip-1) + fnc(ir,it-1,ip-1))
end subroutine interpo_linear_surface_type0
