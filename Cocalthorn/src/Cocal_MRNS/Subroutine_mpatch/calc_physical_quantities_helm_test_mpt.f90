subroutine calc_physical_quantities_helm_test_mpt
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
!
  call calc_scalar_wave_moment
!
  call copy_def_quantities_to_mpt(impt)
!
end subroutine calc_physical_quantities_helm_test_mpt
