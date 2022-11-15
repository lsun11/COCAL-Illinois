module def_vector_x
  use phys_constant, only : long
  implicit none
  real(long), pointer :: vec_xf(:,:,:,:), vec_xg(:,:,:,:)
  real(long), pointer :: hvec_xf(:,:,:,:), hvec_xg(:,:,:,:)
contains
subroutine allocate_vector_x
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use make_array_4d
  implicit none
  call alloc_array4d(vec_xg, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(vec_xf, 0, nrf, 0, ntf, 0, npf, 1, 3)
  call alloc_array4d(hvec_xg, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(hvec_xf, 1, nrf, 1, ntf, 1, npf, 1, 3)
end subroutine allocate_vector_x

subroutine allocate_BHT_vector_x
  use grid_parameter, only : nrg, ntg, npg
  use make_array_4d
  implicit none
  call alloc_array4d(vec_xg, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(hvec_xg, 1, nrg, 1, ntg, 1, npg, 1, 3)
end subroutine allocate_BHT_vector_x
end module def_vector_x
