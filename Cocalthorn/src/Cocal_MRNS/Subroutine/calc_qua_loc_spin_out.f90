subroutine calc_qua_loc_spin_out
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : hrg, rg
  use def_matter_parameter, only : radi
  use def_matter, only : rs
  use make_array_1d
  use make_array_2d
  use def_quantities, only : qua_loc_spin
  use interface_source_qua_loc_spin_out
  use interface_surf_int_grav
  implicit none
  real(long) :: fac8pi, qls
  real(long) :: int_qua_loc_spin
  real(long),pointer :: soug_qua_loc_spin(:,:)
  integer ::  irg, irs, ii
!
  call alloc_array2d(soug_qua_loc_spin, 1, ntg, 1, npg)
  fac8pi = 0.125d0/pi

  do irs=nrf+1, nrf+10
    soug_qua_loc_spin = 0.0d0
    call source_qua_loc_spin_out(soug_qua_loc_spin, irs)
    call surf_int_grav(soug_qua_loc_spin, int_qua_loc_spin, irs)
!
    qls = fac8pi*radi**2*int_qua_loc_spin
    write(6,'(a4,i3,a8,1p,e23.15,a30,1p,e23.15)') 'irs=', irs, '    hrg=', &
    &                                            hrg(irs), 'Quasi local spin NS =', qls

    if (irs==nrf+1)     qua_loc_spin = qls
  end do

  deallocate(soug_qua_loc_spin)

end subroutine calc_qua_loc_spin_out
