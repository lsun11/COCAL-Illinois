subroutine calc_qua_loc_spin_BBH_CF
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg
  use make_array_1d
  use make_array_2d
  use def_quantities, only : qua_loc_spin
  use interface_source_qua_loc_spin
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac8pi
  real(long) :: int_qua_loc_spin
  real(long),pointer :: sou_qua_loc_spin(:,:)
  integer ::  irg
!
  call alloc_array2d(sou_qua_loc_spin, 1, ntg, 1, npg)
!
  fac8pi = 0.125d0/pi

  call source_qua_loc_spin(sou_qua_loc_spin)

  call surf_int_grav_rg(sou_qua_loc_spin, int_qua_loc_spin, 0)
!
  qua_loc_spin = fac8pi*int_qua_loc_spin
!  write(6,'(a28,1p,e14.6)') 'Quasi local spin BH        =', qua_loc_spin

  deallocate(sou_qua_loc_spin)

end subroutine calc_qua_loc_spin_BBH_CF
