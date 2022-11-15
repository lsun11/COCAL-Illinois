subroutine calc_physical_quantities_BBH_trpPunc_CF_mpt
  use phys_constant, only  :  long, nmpt
  use grid_parameter, only : iter_max
  use def_metric, only     : psi
  use def_quantities, only : admmom_asymp
  use interface_calc_fnc_moment_asympto
  implicit none
  integer :: impt, iter_count
  real(long) :: psi_moment(3)
!
! -- Asymptotic mass and angular momentum
  impt = nmpt
  call copy_from_mpatch_all_BBH_CF(impt)
  call copy_def_metric_from_mpt(impt)
  call calc_vector_x_grav(2)
  call calc_vector_phi_grav(2)
  call calc_vector_bh(2)
!
  call excurve_TrpBH_mpt(impt)
  call excurve_TrpBH_gridpoint_mpt(impt)
  call copy_Aij_pBH_to_tfkij
!
  call calc_mass_asympto('bh')
  call calc_angmom_asympto('bh')
  call calc_admmom_asympto('bh')
  call calc_vector_x_grav(1)
  call calc_fnc_moment_asympto(psi,psi_moment)
  admmom_asymp(2) = psi_moment(1)
  call copy_def_quantities_to_mpt(impt)
  call calc_vector_x_grav(2)
!
  do impt = 1, 2
    call copy_from_mpatch_all_BBH_CF(impt)
    call copy_def_metric_from_mpt(impt)
    call calc_vector_x_grav(1)
    call calc_vector_phi_grav(1)
    call calc_vector_bh(1)
    call excurve_TrpBH_mpt(impt)
    call excurve_TrpBH_gridpoint_mpt(impt)
    call copy_Aij_pBH_to_tfkij
    call copy_def_horizon_from_mpt(impt)
    call iteration_AHfinder(iter_count)
    if (iter_count.ge.iter_max) then
      write(6,*)' ** Solution did not converge AH finder **'
    end if
    call calc_AHarea_AHfinder
    call copy_def_quantities_bh_to_mpt(impt)
    call copy_def_horizon_to_mpt(impt)
  end do
!
end subroutine calc_physical_quantities_BBH_trpPunc_CF_mpt
