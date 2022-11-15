! linear interpolation on the GR coordinate.
subroutine interpo_linear1p_type0_2Dsurf(val,fnc)
  use phys_constant, only : long
  implicit none
  real(long) :: val
  real(long) :: fnc(2,2)
  val = 0.25d0*(fnc(1,1) + fnc(1,2)  &
              + fnc(2,1) + fnc(2,2))
end subroutine interpo_linear1p_type0_2Dsurf
