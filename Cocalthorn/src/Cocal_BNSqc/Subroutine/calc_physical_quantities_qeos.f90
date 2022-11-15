subroutine calc_physical_quantities_qeos
  implicit none
  call calc_rest_mass_qeos
  call calc_mass_qeos
  call calc_proper_mass_qeos
  call calc_ang_mom_qeos
  call calc_ToverW_qeos
  call calc_virial_CF_qeos
  call calc_radius
  call calc_redblue_shift
  call calc_mass_asympto('ns')
!  call calc_ang_mom_asymp('ns')
  call excurve_CF_gridpoint
  call calc_angmom_asympto('ns')
  call calc_quad_pole_qeos
  call calc_physq_center_qeos
  call calc_physq_cgs_peos
  call calc_enthalpy_xyzaxis_qeos
  call calc_ergo
end subroutine calc_physical_quantities_qeos
