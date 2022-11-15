!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine copy_weight_midpoint_grav_from_mpt(impt)
  use phys_constant,  only : nnrg, nntg, nnpg
  use grid_parameter, only : nrg, ntg, npg
  use weight_midpoint_grav
  use weight_midpoint_grav_mpt
  use copy_array_static_2dto1d_mpt
  use copy_array_3dto2d_mpt
  use copy_array_4dto3d_mpt
  implicit none
! weight
  integer :: impt
  call copy_arraystatic_2dto1d_mpt(impt, hwdrg_, hwdrg, 1, nnrg)
  call copy_arraystatic_2dto1d_mpt(impt, hwdtg_, hwdtg, 1, nntg)
  call copy_arraystatic_2dto1d_mpt(impt, hwdpg_, hwdpg, 1, nnpg)
  call copy_arraystatic_2dto1d_mpt(impt, tzwdrg_, tzwdrg, 0, nnrg)
  call copy_arraystatic_2dto1d_mpt(impt, tzwdtg_, tzwdtg, 0, nntg)
  call copy_arraystatic_2dto1d_mpt(impt, tzwdpg_, tzwdpg, 0, nnpg)
  call copy_arraystatic_2dto1d_mpt(impt, wdxg_, wdxg, 0, nnrg)
  call copy_array3dto2d_mpt(impt, hwtpgsf_, hwtpgsf, 1, ntg, 1, npg)
  call copy_array4dto3d_mpt(impt, hwrtpg_, hwrtpg, 1, nrg, 1, ntg, 1, npg)
  call copy_array4dto3d_mpt(impt, tzwrtpg_, tzwrtpg, 0, nrg, 1, ntg, 1, npg)
!
end subroutine copy_weight_midpoint_grav_from_mpt
