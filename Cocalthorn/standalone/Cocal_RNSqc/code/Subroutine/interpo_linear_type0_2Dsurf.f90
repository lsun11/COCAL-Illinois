! linear interpolation on the GR coordinate.
subroutine interpo_linear_type0_2Dsurf(val,fnc,it,ip)
  use phys_constant, only : long
  implicit none
  real(long) :: val
  real(long), pointer :: fnc(:,:)
  integer :: it, ip
  val = 0.25d0*(fnc(it,ip)   + fnc(it-1,ip)  &
              + fnc(it,ip-1) + fnc(it-1,ip-1))
end subroutine interpo_linear_type0_2Dsurf
