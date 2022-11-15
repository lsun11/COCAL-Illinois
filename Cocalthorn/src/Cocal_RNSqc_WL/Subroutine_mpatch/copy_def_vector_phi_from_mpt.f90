subroutine copy_def_vector_phi_from_mpt(impt)
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_vector_phi
  use def_vector_phi_mpt
  use copy_array_4dto3d_mpt
  use copy_array_5dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array5dto4d_mpt(impt,vec_phif_,vec_phif, 0,nrf, 0,ntf, 0,npf, 1,3)
  call copy_array5dto4d_mpt(impt,vec_phig_,vec_phig, 0,nrg, 0,ntg, 0,npg, 1,3)
  call copy_array5dto4d_mpt(impt,hvec_phig_,hvec_phig,1,nrg,1,ntg, 1,npg, 1,3)
  call copy_array5dto4d_mpt(impt,hvec_phif_,hvec_phif,1,nrf,1,ntf, 1,npf, 1,3)
  call copy_array4dto3d_mpt(impt,hvec_phif_surface_, hvec_phif_surface, &
  &                                                         1,ntf, 1,npf, 1,3)
!
end subroutine copy_def_vector_phi_from_mpt
