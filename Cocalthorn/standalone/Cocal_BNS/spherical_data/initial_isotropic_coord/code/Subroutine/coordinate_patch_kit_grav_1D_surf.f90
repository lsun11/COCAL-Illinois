subroutine coordinate_patch_kit_grav_1D_surf(igrid)
  use grid_parameter_1D
  use coordinate_grav_r_1D
  use weight_grav_1D
  use radial_perm_fn_grav_1D
  implicit none
  integer :: igrid
! call subroutines. the order is important.
  call read_parameter_1D
  call read_parameter_1D_surf

  if (igrid==3) then
    write(6,*) "***Neutron star grid with rg(nrf)<1... ****"
    call grid_r_1D_surf
  else if (igrid==4) then
    write(6,*) "***Neutron star grid with rg(nrf)=1 and constant dr until r=7...****"
    call grid_r_1D_surf_const
  else
    write(6,*) "***Choose a grid: 3 or 4...exiting****"
    stop
  end if

  call calc_hfsn
  call weight_calc_grav_1D
end subroutine coordinate_patch_kit_grav_1D_surf
