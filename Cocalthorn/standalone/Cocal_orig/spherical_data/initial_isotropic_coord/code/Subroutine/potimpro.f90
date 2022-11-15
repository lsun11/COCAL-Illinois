subroutine potimpro(potg,backg,fffac)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrg
  implicit none
!
  real(8), intent(inout) :: potg(0:nnrg), backg(0:nnrg), fffac
!
  potg(0:nrg) = (1.0d0 - fffac)*backg(0:nrg) + fffac * potg(0:nrg)
!
end subroutine potimpro
