subroutine copy_def_vector_phi_to_mpt(impt)
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_vector_phi
  use def_vector_phi_mpt
  use copy_array_3dto4d_mpt
  use copy_array_4dto5d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto5d_mpt(impt, vec_phig, vec_phig_, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call copy_array4dto5d_mpt(impt, vec_phif, vec_phif_, 0, nrf, 0, ntf, 0, npf, 1, 3)
  call copy_array4dto5d_mpt(impt, hvec_phig, hvec_phig_, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call copy_array4dto5d_mpt(impt, hvec_phif, hvec_phif_, 1, nrf, 1, ntf, 1, npf, 1, 3)
  call copy_array3dto4d_mpt(impt, hvec_phif_surface, hvec_phif_surface_, 1, ntf, 1, npf, 1, 3)
!
end subroutine copy_def_vector_phi_to_mpt
