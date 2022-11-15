subroutine calc_physical_quantities_BBH_CF_mpt
  use phys_constant, only :  long, nmpt
  implicit none
  integer :: impt
!
! -- Asymptotic mass and angular momentum
  impt = nmpt
  call copy_from_mpatch_all_BBH_CF(impt)
  call copy_def_metric_from_mpt(impt)
  call calc_vector_x_grav(2)
  call calc_vector_phi_grav(2)
  call calc_vector_bh(2)
  call excurve_CF_gridpoint_bhex
  call calc_mass_asympto('bh')
  call calc_angmom_asympto('bh')
  call calc_admmom_asympto('bh')
  call copy_def_quantities_to_mpt(impt)
!
  do impt = 1, 2
    call copy_from_mpatch_all_BBH_CF(impt)
    call copy_def_metric_from_mpt(impt)
    call calc_AH_BBH_CF
    call copy_def_quantities_bh_to_mpt(impt)
  end do
!
end subroutine calc_physical_quantities_BBH_CF_mpt
