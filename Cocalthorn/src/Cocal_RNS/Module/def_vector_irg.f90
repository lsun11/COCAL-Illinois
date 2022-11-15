module def_vector_irg
  use phys_constant, only : long
  implicit none
  real(long), pointer ::  vec_irg_cm_xg(:,:,:),  vec_irg_cm_phig(:,:,:)
  real(long), pointer :: hvec_irg_cm_xg(:,:,:), hvec_irg_cm_phig(:,:,:)
  real(long), pointer ::  vec_irg_cbh_xg(:,:,:), vec_irg_cbh_phig(:,:,:)
  real(long), pointer :: hvec_irg_cbh_xg(:,:,:),hvec_irg_cbh_phig(:,:,:)
!  real(long), pointer ::  vec_irg_cbh_spin(:,:,:)
!  real(long), pointer :: hvec_irg_cbh_spin(:,:,:)
contains
subroutine allocate_vector_irg
  use grid_parameter, only : ntg, npg
  use make_array_3d
  implicit none
  call alloc_array3d( vec_irg_cm_xg, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_irg_cm_xg, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_irg_cm_phig, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_irg_cm_phig, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_irg_cbh_xg, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_irg_cbh_xg, 1, ntg, 1, npg, 1, 3)
  call alloc_array3d( vec_irg_cbh_phig, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(hvec_irg_cbh_phig, 1, ntg, 1, npg, 1, 3)
!  call alloc_array3d( vec_irg_cbh_spin, 0, ntg, 0, npg, 1, 3)
!  call alloc_array3d(hvec_irg_cbh_spin, 1, ntg, 1, npg, 1, 3)
end subroutine allocate_vector_irg
end module def_vector_irg
