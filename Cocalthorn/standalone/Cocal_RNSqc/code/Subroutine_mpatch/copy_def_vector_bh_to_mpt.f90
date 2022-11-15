subroutine copy_def_vector_bh_to_mpt(impt)
  use grid_parameter, only : ntg, npg
  use def_vector_bh
  use def_vector_bh_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt, vec_bh_cm_xg,    vec_bh_cm_xg_,    0,ntg,0,npg,1,3)
  call copy_array3dto4d_mpt(impt,hvec_bh_cm_xg,   hvec_bh_cm_xg_,    1,ntg,1,npg,1,3)
  call copy_array3dto4d_mpt(impt, vec_bh_cm_phig,  vec_bh_cm_phig_,  0,ntg,0,npg,1,3)
  call copy_array3dto4d_mpt(impt,hvec_bh_cm_phig, hvec_bh_cm_phig_,  1,ntg,1,npg,1,3)
  call copy_array3dto4d_mpt(impt, vec_bh_cbh_xg,   vec_bh_cbh_xg_,   0,ntg,0,npg,1,3)
  call copy_array3dto4d_mpt(impt,hvec_bh_cbh_xg,  hvec_bh_cbh_xg_,   1,ntg,1,npg,1,3)
  call copy_array3dto4d_mpt(impt, vec_bh_cbh_phig, vec_bh_cbh_phig_, 0,ntg,0,npg,1,3)
  call copy_array3dto4d_mpt(impt,hvec_bh_cbh_phig,hvec_bh_cbh_phig_, 1,ntg,1,npg,1,3)
  call copy_array3dto4d_mpt(impt, vec_bh_cbh_spin, vec_bh_cbh_spin_, 0,ntg,0,npg,1,3)
  call copy_array3dto4d_mpt(impt,hvec_bh_cbh_spin,hvec_bh_cbh_spin_, 1,ntg,1,npg,1,3)
!
end subroutine copy_def_vector_bh_to_mpt
