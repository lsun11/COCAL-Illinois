! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
subroutine copy_radial_green_fn_grav_from_mpt(impt)
  use grid_parameter, only : nrg, nlg
  use radial_green_fn_grav
  use radial_green_fn_grav_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, hgfn_, hgfn, 1, nrg, 0, nlg, 0, nrg)
  call copy_array4dto3d_mpt(impt, gfnsf_, gfnsf, 0, nlg, 0, nrg, 1, 4)
!
end subroutine copy_radial_green_fn_grav_from_mpt
