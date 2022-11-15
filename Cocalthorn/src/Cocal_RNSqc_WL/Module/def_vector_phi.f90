module def_vector_phi
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  implicit none
  real(long), pointer :: vec_phif(:,:,:,:), vec_phig(:,:,:,:)
  real(long), pointer :: hvec_phif(:,:,:,:), hvec_phig(:,:,:,:)
  real(long), pointer :: hvec_phif_surface(:,:,:)
contains
subroutine allocate_vector_phi
  use make_array_3d
  use make_array_4d
  implicit none
  call alloc_array4d(vec_phig, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(vec_phif, 0, nrf, 0, ntf, 0, npf, 1, 3)
  call alloc_array4d(hvec_phig, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(hvec_phif, 1, nrf, 1, ntf, 1, npf, 1, 3)
  call alloc_array3d(hvec_phif_surface, 1, ntf, 1, npf, 1, 3)
end subroutine allocate_vector_phi

subroutine allocate_BHT_vector_phi
  use make_array_3d
  use make_array_4d
  implicit none
  call alloc_array4d(vec_phig, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(hvec_phig, 1, nrg, 1, ntg, 1, npg, 1, 3)
end subroutine allocate_BHT_vector_phi
end module def_vector_phi
