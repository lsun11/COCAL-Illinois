subroutine copy_def_vector_bh_from_mpt(impt)
  use grid_parameter, only : ntg, npg
  use def_vector_bh
  use def_vector_bh_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, vec_bh_cm_xg_,    vec_bh_cm_xg,    0,ntg,0,npg,1,3)
  call copy_array4dto3d_mpt(impt,hvec_bh_cm_xg_,   hvec_bh_cm_xg,    1,ntg,1,npg,1,3)
  call copy_array4dto3d_mpt(impt, vec_bh_cm_phig_,  vec_bh_cm_phig,  0,ntg,0,npg,1,3)
  call copy_array4dto3d_mpt(impt,hvec_bh_cm_phig_, hvec_bh_cm_phig,  1,ntg,1,npg,1,3)
  call copy_array4dto3d_mpt(impt, vec_bh_cbh_xg_,   vec_bh_cbh_xg,   0,ntg,0,npg,1,3)
  call copy_array4dto3d_mpt(impt,hvec_bh_cbh_xg_,  hvec_bh_cbh_xg,   1,ntg,1,npg,1,3)
  call copy_array4dto3d_mpt(impt, vec_bh_cbh_phig_, vec_bh_cbh_phig, 0,ntg,0,npg,1,3)
  call copy_array4dto3d_mpt(impt,hvec_bh_cbh_phig_,hvec_bh_cbh_phig, 1,ntg,1,npg,1,3)
  call copy_array4dto3d_mpt(impt, vec_bh_cbh_spin_, vec_bh_cbh_spin, 0,ntg,0,npg,1,3)
  call copy_array4dto3d_mpt(impt,hvec_bh_cbh_spin_,hvec_bh_cbh_spin, 1,ntg,1,npg,1,3)
!
end subroutine copy_def_vector_bh_from_mpt
