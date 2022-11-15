subroutine gr2fl(grv,flv)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrf
  implicit none
  real(8), intent(inout) :: flv(0:nnrg), grv(0:nnrg)
!
  flv(0:nrf) = grv(0:nrf)
!
end subroutine gr2fl
