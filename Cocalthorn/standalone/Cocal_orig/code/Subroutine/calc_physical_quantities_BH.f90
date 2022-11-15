subroutine calc_physical_quantities_BH
  implicit none
!
  call calc_vector_x_grav(1)
  call calc_vector_phi_grav(1)
  call copy_Aij_pBH_to_tfkij
  call calc_mass_pBH
  call calc_mass_asympto('bh')
  call calc_angmom_asympto('bh')
  call calc_admmom_asympto('bh')
!
end subroutine calc_physical_quantities_BH
