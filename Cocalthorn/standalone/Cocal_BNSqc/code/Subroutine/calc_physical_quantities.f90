subroutine calc_physical_quantities
  implicit none
  call calc_rest_mass
  call calc_mass
  call calc_proper_mass
  call calc_ang_mom
  call calc_ToverW
  call calc_radius
  call calc_redblue_shift
  call calc_mass_asympto('ns')
  call calc_ang_mom_asymp('ns')
  call calc_physq_center
end subroutine calc_physical_quantities
