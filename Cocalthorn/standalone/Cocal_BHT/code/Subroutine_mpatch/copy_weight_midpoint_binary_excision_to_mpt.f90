!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine copy_weight_midpoint_binary_excision_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use weight_midpoint_binary_excision
  use weight_midpoint_binary_excision_mpt
  use copy_array_2dto3d_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
! weight binary excision
  call copy_array3dto4d_mpt(impt,hwrtpg_ex,hwrtpg_ex_,1,nrg,1,ntg,1,npg)
  call copy_array2dto3d_mpt(impt,hwtpg_ex, hwtpg_ex_ ,1,ntg,1,npg)
!
end subroutine copy_weight_midpoint_binary_excision_to_mpt
