module def_vector_bh_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer ::  vec_bh_cm_xg_(:,:,:,:),  vec_bh_cm_phig_(:,:,:,:)
  real(long), pointer :: hvec_bh_cm_xg_(:,:,:,:), hvec_bh_cm_phig_(:,:,:,:)
  real(long), pointer ::  vec_bh_cbh_xg_(:,:,:,:), vec_bh_cbh_phig_(:,:,:,:)
  real(long), pointer :: hvec_bh_cbh_xg_(:,:,:,:),hvec_bh_cbh_phig_(:,:,:,:)
  real(long), pointer ::  vec_bh_cbh_spin_(:,:,:,:)
  real(long), pointer :: hvec_bh_cbh_spin_(:,:,:,:)
end module def_vector_bh_mpt
