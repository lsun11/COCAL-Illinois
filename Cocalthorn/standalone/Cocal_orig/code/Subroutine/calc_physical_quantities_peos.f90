subroutine calc_physical_quantities_peos
  implicit none
  call calc_rest_mass_peos
  call calc_mass_peos
  call calc_proper_mass_peos
  call calc_ang_mom_peos
  call calc_ToverW_peos
  call calc_virial_CF
  call calc_radius
  call calc_redblue_shift
  call calc_mass_asympto('ns')
!  call calc_ang_mom_asymp('ns')
  call excurve_CF_gridpoint
  call calc_angmom_asympto('ns')
  call calc_quad_pole_peos
  call calc_physq_center_peos
  call calc_physq_cgs_peos
  call calc_enthalpy_xyzaxis
end subroutine calc_physical_quantities_peos
