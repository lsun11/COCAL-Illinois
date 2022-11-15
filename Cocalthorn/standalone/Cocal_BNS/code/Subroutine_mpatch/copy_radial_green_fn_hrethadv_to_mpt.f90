subroutine copy_radial_green_fn_hrethadv_to_mpt(impt)
  use radial_green_fn_hrethadv
  use radial_green_fn_hrethadv_mpt
  use grid_parameter, only : nrg, nlg
  use copy_array_3dto4d_mpt
  use copy_array_4dto5d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto5d_mpt(impt, bsjy, bsjy_, 1, nrg, 0, nlg, 0, nlg, 0, nrg)
  call copy_array3dto4d_mpt(impt, sbsjy, sbsjy_, 0, nlg, 0, nlg, 0, nrg)
  call copy_array3dto4d_mpt(impt, sbsjyp, sbsjyp_, 0, nlg, 0, nlg, 0, nrg)
!
end subroutine copy_radial_green_fn_hrethadv_to_mpt
