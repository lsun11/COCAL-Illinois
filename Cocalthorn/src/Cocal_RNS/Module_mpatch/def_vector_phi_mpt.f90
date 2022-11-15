module def_vector_phi_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer :: vec_phif_(:,:,:,:,:), vec_phig_(:,:,:,:,:)
  real(long), pointer :: hvec_phif_(:,:,:,:,:), hvec_phig_(:,:,:,:,:)
  real(long), pointer :: hvec_phif_surface_(:,:,:,:)
end module def_vector_phi_mpt
