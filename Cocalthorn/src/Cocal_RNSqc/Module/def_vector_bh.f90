module def_vector_bh
  use phys_constant, only : long
  implicit none
  real(long), pointer ::  vec_bh_cm_xg(:,:,:),  vec_bh_cm_phig(:,:,:)
  real(long), pointer :: hvec_bh_cm_xg(:,:,:), hvec_bh_cm_phig(:,:,:)
  real(long), pointer ::  vec_bh_cbh_xg(:,:,:), vec_bh_cbh_phig(:,:,:)
  real(long), pointer :: hvec_bh_cbh_xg(:,:,:),hvec_bh_cbh_phig(:,:,:)
  real(long), pointer ::  vec_bh_cbh_spin(:,:,:)
  real(long), pointer :: hvec_bh_cbh_spin(:,:,:)
contains
subroutine allocate_vector_bh
  use grid_parameter, only : ntg, npg
  use make_array_3d
  implicit none
  call alloc_array3d( vec_bh_cm_xg, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_bh_cm_xg, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_bh_cm_phig, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_bh_cm_phig, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_bh_cbh_xg, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_bh_cbh_xg, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_bh_cbh_phig, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_bh_cbh_phig, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_bh_cbh_spin, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_bh_cbh_spin, 1, ntg, 1, npg, 1, 3)
end subroutine allocate_vector_bh
end module def_vector_bh
