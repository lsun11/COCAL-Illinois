subroutine allocate_excurve_on_SFC_CF
  use grid_parameter, only : nrf, ntf, npf
  use def_excurve_on_SFC_CF
  use make_array_5d
  implicit none
!
  call alloc_array5d(tfkij_grid_fluid, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, 3)
  tfkij_grid_fluid(0:nrf,0:ntf,0:npf,1:3,1:3) = 0.0d0
!
end subroutine allocate_excurve_on_SFC_CF
