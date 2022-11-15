!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine copy_weight_midpoint_grav_to_mpt(impt)
  use phys_constant,  only : nnrg, nntg, nnpg
  use grid_parameter, only : nrg, ntg, npg
  use weight_midpoint_grav
  use weight_midpoint_grav_mpt
  use copy_array_static_1dto2d_mpt
  use copy_array_2dto3d_mpt
  use copy_array_3dto4d_mpt
  implicit none
! weight
  integer :: impt
  call copy_arraystatic_1dto2d_mpt(impt, hwdrg, hwdrg_, 1, nnrg)
  call copy_arraystatic_1dto2d_mpt(impt, hwdtg, hwdtg_, 1, nntg)
  call copy_arraystatic_1dto2d_mpt(impt, hwdpg, hwdpg_, 1, nnpg)
  call copy_arraystatic_1dto2d_mpt(impt, tzwdrg, tzwdrg_, 0, nnrg)
  call copy_arraystatic_1dto2d_mpt(impt, tzwdtg, tzwdtg_, 0, nntg)
  call copy_arraystatic_1dto2d_mpt(impt, tzwdpg, tzwdpg_, 0, nnpg)
  call copy_arraystatic_1dto2d_mpt(impt, wdxg, wdxg_, 0, nnrg)
  call copy_array2dto3d_mpt(impt, hwtpgsf, hwtpgsf_, 1, ntg, 1, npg)
  call copy_array3dto4d_mpt(impt, hwrtpg, hwrtpg_, 1, nrg, 1, ntg, 1, npg)
  call copy_array3dto4d_mpt(impt, tzwrtpg, tzwrtpg_, 0, nrg, 1, ntg, 1, npg)
!
end subroutine copy_weight_midpoint_grav_to_mpt
