subroutine next_solution
  use phys_constant,  only : long
  use grid_parameter, only : nrf_deform, deform_par, nrf, ratio
  implicit none
  nrf_deform = nrf_deform + deform_par
  ratio = dble(nrf_deform)/dble(nrf)
end subroutine next_solution
