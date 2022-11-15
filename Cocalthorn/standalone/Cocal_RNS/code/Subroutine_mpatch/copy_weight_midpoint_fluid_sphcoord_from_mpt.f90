!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine copy_weight_midpoint_fluid_sphcoord_from_mpt(impt)
  use grid_parameter, only : nrf, ntg, npg
  use weight_midpoint_fluid_sphcoord
  use weight_midpoint_fluid_sphcoord_mpt
  implicit none
  integer :: impt
! weight binary excision
  call copy_array4dto3d_mpt(impt,hwrtpg_fc_,hwrtpg_fc,1,nrf,1,ntg,1,npg)
  call copy_array3dto2d_mpt(impt, hwtpg_fc_, hwtpg_fc,1,ntg,1,npg)
!
end subroutine copy_weight_midpoint_fluid_sphcoord_from_mpt
