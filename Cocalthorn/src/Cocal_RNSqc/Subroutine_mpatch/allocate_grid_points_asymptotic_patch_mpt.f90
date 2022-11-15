subroutine allocate_grid_points_asymptotic_patch_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_asymptotic_patch_mpt
  use make_array_4d
  implicit none
!
! -- ra, tha, phia, are extended coordinates ONLY in radial direction.
!
  call alloc_array4d(ra_, -2, nrg+2, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(tha_, -2, nrg+2, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(phia_, -2, nrg+2, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hra_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array4d(htha_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array4d(hphia_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
!
end subroutine allocate_grid_points_asymptotic_patch_mpt
