subroutine allocate_grid_points_binary_in_asympto_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_binary_in_asympto_mpt
  use make_array_4d
  implicit none
!
! -- ra, tha, phia, are extended coordinates ONLY in radial direction.
!
  call alloc_array4d(rb_a_, -2, nrg+2, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(thb_a_, -2, nrg+2, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(phib_a_, -2, nrg+2, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hrb_a_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array4d(hthb_a_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array4d(hphib_a_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
!
end subroutine allocate_grid_points_binary_in_asympto_mpt
