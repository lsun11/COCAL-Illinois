module def_vector_x_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer :: vec_xf_(:,:,:,:,:), vec_xg_(:,:,:,:,:)
  real(long), pointer :: hvec_xf_(:,:,:,:,:), hvec_xg_(:,:,:,:,:)
end module def_vector_x_mpt
