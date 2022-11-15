!  Associated Legendre function and factorials
!______________________________________________
subroutine allocate_legendre_fn_irbns_grav_mpt
  use legendre_fn_irbns_grav_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : ntg, nlg
  use make_array_2d
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array4d(plmg_, 0, nlg, 0, nlg, 0, ntg, 1, nmpt)
  call alloc_array4d(hplmg_, 0, nlg, 0, nlg, 1, ntg, 1, nmpt)
  call alloc_array4d(dtplmg_, 0, nlg, 0, nlg, 0, ntg, 1, nmpt)
  call alloc_array4d(hdtplmg_, 0, nlg, 0, nlg, 1, ntg, 1, nmpt)
  call alloc_array4d(yplmg_, 0, nlg, 0, nlg, 0, ntg, 1, nmpt)
  call alloc_array4d(hyplmg_, 0, nlg, 0, nlg, 1, ntg, 1, nmpt)
  call alloc_array4d(dtyplmg_, 0, nlg, 0, nlg, 0, ntg, 1, nmpt)
  call alloc_array4d(hdtyplmg_, 0, nlg, 0, nlg, 1, ntg, 1, nmpt)
  call alloc_array3d(facnmg_, 0, nlg, 0, nlg, 1, nmpt)
  call alloc_array2d(epsig_, 0, nlg, 1, nmpt)
!
end subroutine allocate_legendre_fn_irbns_grav_mpt
