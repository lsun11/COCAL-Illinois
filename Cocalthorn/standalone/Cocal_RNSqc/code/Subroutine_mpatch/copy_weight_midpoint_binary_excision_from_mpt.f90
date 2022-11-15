!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine copy_weight_midpoint_binary_excision_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use weight_midpoint_binary_excision
  use weight_midpoint_binary_excision_mpt
  use copy_array_3dto2d_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
! weight binary excision
  call copy_array4dto3d_mpt(impt,hwrtpg_ex_,hwrtpg_ex,1,nrg,1,ntg,1,npg)
  call copy_array3dto2d_mpt(impt,hwtpg_ex_, hwtpg_ex, 1,ntg,1,npg)
!
end subroutine copy_weight_midpoint_binary_excision_from_mpt
