subroutine copy_def_quantities_derived_to_mpt(impt)
  use def_quantities_derived
  use def_quantities_derived_mpt
  
  implicit none
  integer :: i, impt
!
  i=0
  i=i+1; def_quantities_derived_real_(i,impt) = coord_radius_normalized_x
  i=i+1; def_quantities_derived_real_(i,impt) = coord_radius_normalized_y
  i=i+1; def_quantities_derived_real_(i,impt) = coord_radius_normalized_z
  i=i+1; def_quantities_derived_real_(i,impt) = coord_axis_ratio_yx
  i=i+1; def_quantities_derived_real_(i,impt) = proper_axis_ratio_yx
  i=i+1; def_quantities_derived_real_(i,impt) = coord_axis_ratio_zx
  i=i+1; def_quantities_derived_real_(i,impt) = proper_axis_ratio_zx
  i=i+1; def_quantities_derived_real_(i,impt) = coord_eccentricity_zx
  i=i+1; def_quantities_derived_real_(i,impt) = proper_eccentricity_zx
  i=i+1; def_quantities_derived_real_(i,impt) = omega
  i=i+1; def_quantities_derived_real_(i,impt) = omega_M
  i=i+1; def_quantities_derived_real_(i,impt) = omega2_M2
  i=i+1; def_quantities_derived_real_(i,impt) = E_bind_over_M
  i=i+1; def_quantities_derived_real_(i,impt) = J_over_M2
  i=i+1; def_quantities_derived_real_(i,impt) = J_over_Madm2
  i=i+1; def_quantities_derived_real_(i,impt) = Madm_MK_virial
  i=i+1; def_quantities_derived_real_(i,impt) = M0_minus_M0sph_normalized
  i=i+1; def_quantities_derived_real_(i,impt) = rho_c_normalized
  i=i+1; def_quantities_derived_real_(i,impt) = pre_c_normalized
  i=i+1; def_quantities_derived_real_(i,impt) = epsilon_c_normalized
  i=i+1; def_quantities_derived_real_(i,impt) = p_over_rho_c_normalized

  
!
end subroutine copy_def_quantities_derived_to_mpt
