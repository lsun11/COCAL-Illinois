subroutine fl2gr(flv,grv)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrf, nrg
  implicit none
  real(8), intent(inout) :: flv(0:nnrg), grv(0:nnrg)
!
  grv(0:nrg) = 0.0d0
  grv(0:nrf) = flv(0:nrf)
!
end subroutine fl2gr
