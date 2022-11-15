subroutine calc_qua_loc_spin_peos
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : hrg, rg
  use def_matter_parameter, only : radi
  use def_matter, only : rs
  use make_array_1d
  use make_array_2d
  use def_quantities, only : qua_loc_spin, qua_loc_spin_surf
  use interface_source_qua_loc_spin_peos_fluid
  use interface_source_qua_loc_spin_peos
  use interface_surf_int_grav
  use interface_surf_int_fluid_rg
  implicit none
  real(long) :: fac8pi, qls
  real(long) :: int_qua_loc_spin
  real(long),pointer :: soug_qua_loc_spin(:,:), souf_qua_loc_spin(:,:)
  integer ::  irg, irs, ii
!
  call alloc_array2d(soug_qua_loc_spin, 1, ntg, 1, npg)
  call alloc_array2d(souf_qua_loc_spin, 1, ntf, 1, npf)
!
  fac8pi = 0.125d0/pi

  do irs = nrf-5, nrf
    call source_qua_loc_spin_peos_fluid(souf_qua_loc_spin, irs)
    call surf_int_fluid_rg(souf_qua_loc_spin, int_qua_loc_spin, irs)
!
    qls = fac8pi*radi**2*int_qua_loc_spin
    write(6,'(a4,i3,a8,1p,e23.15,a30,1p,e23.15)') 'irs=', irs, '  rs*rg=', &
    &                               rs(ntf/2,0)*rg(irs), 'Quasi local spin NS =', qls
    
    if (irs==nrf)  qua_loc_spin_surf = qls
  end do
 
  call calc_vector_x_grav(1)
  call calc_vector_phi_grav(1)

  do irs=nrf+1, nrf+10
    soug_qua_loc_spin = 0.0d0
    call source_qua_loc_spin_peos(soug_qua_loc_spin, irs)
    call surf_int_grav(soug_qua_loc_spin, int_qua_loc_spin, irs)
!
    qls = fac8pi*radi**2*int_qua_loc_spin
    write(6,'(a4,i3,a8,1p,e23.15,a30,1p,e23.15)') 'irs=', irs, '    hrg=', &
    &                                            hrg(irs), 'Quasi local spin NS =', qls

    if (irs==nrf+1)     qua_loc_spin = qls
  end do

  deallocate(soug_qua_loc_spin)
  deallocate(souf_qua_loc_spin)

end subroutine calc_qua_loc_spin_peos
