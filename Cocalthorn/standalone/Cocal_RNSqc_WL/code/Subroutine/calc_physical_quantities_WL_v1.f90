subroutine calc_physical_quantities_WL_v1
  implicit none
!
  call interpo_gr2fl_metric_WL
!
  call calc_rest_mass_peos
  call calc_mass_WL
  call calc_proper_mass_peos
  call calc_ang_mom_WL
  call calc_ToverW_WL   ! to compute T_kinene_omJ, I_inertia etc
  call calc_virial_WL
  call calc_radius_WL
  call calc_redblue_shift_WL
!
  call excurve_WL_gridpoint
  call calc_mass_asympto('ns')
  call calc_angmom_asympto('ns')
!
  call calc_quad_pole_peos
  call calc_physq_center_peos
  call calc_physq_cgs_peos
  call calc_enthalpy_xyzaxis
!
  call calc_circ_line_peos
  call calc_circ_surf_peos
  call calc_soundspeed_peos
end subroutine calc_physical_quantities_WL_v1
