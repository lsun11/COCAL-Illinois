subroutine coordinate_patch_kit_grav_grid(igrid)
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use weight_midpoint_grav
!  use weight_midpoint_fluid
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use radial_green_fn_grav
  use legendre_fn_grav
  implicit none
  integer :: igrid

  if (igrid==1) then
    write(6,*) "Default grid..."
    call grid_r

  else if (igrid==2) then
    write(6,*) "Black hole grid..."
    call grid_r_bhex('eBH')

  else if (igrid==3) then
    write(6,*) "Neutron star grid with rg(nrf)<1..."
    call grid_r_bns

  else if (igrid==4) then
    write(6,*) "Neutron star grid with rg(nrf)=1 and constant dr until r=7..."
    call grid_r_bns_const

  else if (igrid==5) then
    write(6,*) "*** Black hole torus grid... ****"
    call grid_r_bht('eBH')

  else if (igrid==6) then
    write(6,*) "*** Neutron star torus grid... ****"
!    call grid_r_nst('eBH')

  else
    write(6,*) "Choose a grid: 1-6...exiting"
    stop
  end if

!  call grid_r
!  call allocate_hgfn
!  call calc_hgfn
  call grid_theta
  call trig_grav_theta
  call allocate_legendre
  call legendre
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call allocate_weight_midpoint_grav
  call weight_calc_midpoint_grav
!  call weight_calc_midpoint_grav_th4th
!  call allocate_weight_midpoint_fluid
!  call weight_calc_midpoint_fluid
  call grid_extended
end subroutine coordinate_patch_kit_grav_grid

