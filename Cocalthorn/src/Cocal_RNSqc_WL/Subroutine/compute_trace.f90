subroutine compute_trace(met,fnc,val)
  use phys_constant, only : long
  implicit none
  real(long), intent(out) :: val
  real(long), intent(in)  :: met(3,3), fnc(3,3)
  val = met(1,1)*fnc(1,1) + met(1,2)*fnc(1,2) + met(1,3)*fnc(1,3) &
  &   + met(2,1)*fnc(2,1) + met(2,2)*fnc(2,2) + met(2,3)*fnc(2,3) &
  &   + met(3,1)*fnc(3,1) + met(3,2)*fnc(3,2) + met(3,3)*fnc(3,3)
end subroutine compute_trace
