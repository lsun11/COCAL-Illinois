! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
subroutine allocate_radial_green_fn_grav_mpt
  use phys_constant, only : nmpt
  use make_array_4d
  use grid_parameter, only : nrg, nlg
  use radial_green_fn_grav_mpt
  implicit none
!
  call alloc_array4d(hgfn_, 1, nrg, 0, nlg, 0, nrg, 1, nmpt)
  call alloc_array4d(gfnsf_, 0, nlg, 0, nrg, 1, 4, 1, nmpt)
!
end subroutine allocate_radial_green_fn_grav_mpt
