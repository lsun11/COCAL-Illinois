subroutine copy_radial_green_fn_hrethadv_from_mpt(impt)
  use radial_green_fn_hrethadv
  use radial_green_fn_hrethadv_mpt
  use grid_parameter, only : nrg, nlg
  implicit none
  integer :: impt
!
  call copy_array5dto4d_mpt(impt, bsjy_, bsjy, 1, nrg, 0, nlg, 0, nlg, 0, nrg)
  call copy_array4dto3d_mpt(impt, sbsjy_, sbsjy, 0, nlg, 0, nlg, 0, nrg)
  call copy_array4dto3d_mpt(impt, sbsjyp_, sbsjyp, 0, nlg, 0, nlg, 0, nrg)
!
end subroutine copy_radial_green_fn_hrethadv_from_mpt
