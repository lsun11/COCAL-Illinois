subroutine calc_physical_quantities_BBH_CF
  implicit none
  call calc_mass_BBH_CF
  call calc_mass_BBH_CF_adm_vol
  call calc_ang_mom_BBH_CF_inf
!  call modify_r0_excurve
  call calc_ang_mom_BBH_CF_thr
  call calc_ang_mom_BBH_CF_smarr
  call calc_app_hor_area_BBH_CF
  call calc_qua_loc_spin_BBH_CF
end subroutine calc_physical_quantities_BBH_CF
