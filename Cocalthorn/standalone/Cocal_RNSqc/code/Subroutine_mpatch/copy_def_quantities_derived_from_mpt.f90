subroutine copy_def_quantities_derived_from_mpt(impt)
  use def_quantities_derived
  use def_quantities_derived_mpt
  
  implicit none
  integer :: i, impt
!
  i=0
  i=i+1;   coord_radius_normalized_x = def_quantities_derived_real_(i,impt)
  i=i+1;   coord_radius_normalized_y = def_quantities_derived_real_(i,impt)
  i=i+1;   coord_radius_normalized_z = def_quantities_derived_real_(i,impt)
  i=i+1;   coord_axis_ratio_yx    = def_quantities_derived_real_(i,impt)
  i=i+1;   proper_axis_ratio_yx   = def_quantities_derived_real_(i,impt)
  i=i+1;   coord_axis_ratio_zx    = def_quantities_derived_real_(i,impt)
  i=i+1;   proper_axis_ratio_zx   = def_quantities_derived_real_(i,impt)
  i=i+1;   coord_eccentricity_zx  = def_quantities_derived_real_(i,impt)
  i=i+1;   proper_eccentricity_zx = def_quantities_derived_real_(i,impt)
  i=i+1;   omega      = def_quantities_derived_real_(i,impt)
  i=i+1;   omega_M    = def_quantities_derived_real_(i,impt)
  i=i+1;   omega2_M2  = def_quantities_derived_real_(i,impt)
  i=i+1;   E_bind_over_M  = def_quantities_derived_real_(i,impt)
  i=i+1;   J_over_M2      = def_quantities_derived_real_(i,impt)
  i=i+1;   J_over_Madm2   = def_quantities_derived_real_(i,impt)
  i=i+1;   Madm_MK_virial = def_quantities_derived_real_(i,impt)
  i=i+1;   M0_minus_M0sph_normalized = def_quantities_derived_real_(i,impt)
  i=i+1;   rho_c_normalized          = def_quantities_derived_real_(i,impt)
  i=i+1;   pre_c_normalized          = def_quantities_derived_real_(i,impt)
  i=i+1;   epsilon_c_normalized      = def_quantities_derived_real_(i,impt)
  i=i+1;   p_over_rho_c_normalized   = def_quantities_derived_real_(i,impt)

  
!
end subroutine copy_def_quantities_derived_from_mpt
