subroutine allocate_grid_points_binary_excision_mpt
  use grid_points_binary_excision_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use make_array_3d
  use make_array_4d
  use make_int_array_3d
  implicit none
!
  call alloc_int_array3d(irg_exin_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_int_array3d(irg_exout_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_int_array3d(itg_exin_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_int_array3d(itg_exout_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_int_array3d(ipg_exin_, 0, nrg, 0, ntg, 1, nmpt)
  call alloc_int_array3d(ipg_exout_,0, nrg, 0, ntg, 1, nmpt)
  call alloc_array3d(rg_exin_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array3d(rg_exout_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array3d(thg_exin_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_array3d(thg_exout_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_array3d(phig_exin_, 0, nrg, 0, ntg, 1, nmpt)
  call alloc_array3d(phig_exout_, 0, nrg, 0, ntg, 1, nmpt)
!
  call alloc_int_array3d(ihrg_exin_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_int_array3d(ihrg_exout_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_int_array3d(ihtg_exin_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_int_array3d(ihtg_exout_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_int_array3d(ihpg_exin_, 0, nrg, 0, ntg, 1, nmpt)
  call alloc_int_array3d(ihpg_exout_, 0, nrg, 0, ntg, 1, nmpt)
  call alloc_array3d(hrg_exin_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array3d(hrg_exout_, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array3d(hthg_exin_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_array3d(hthg_exout_, 0, nrg, 0, npg, 1, nmpt)
  call alloc_array3d(hphig_exin_, 0, nrg, 0, ntg, 1, nmpt)
  call alloc_array3d(hphig_exout_, 0, nrg, 0, ntg, 1, nmpt)
!
  call alloc_array4d(rb_, 0,nrg,0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(thb_, 0,nrg,0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(phib_, 0,nrg,0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hrb_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array4d(hthb_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array4d(hphib_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
!
end subroutine allocate_grid_points_binary_excision_mpt
