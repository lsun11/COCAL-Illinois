!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_fluid_sphcoord_mpt
  use phys_constant, only : long
  implicit none
! weight binary excision
  real(long), pointer :: hwrtpg_fc_(:,:,:,:)
  real(long), pointer :: hwtpg_fc_(:,:,:)
end module weight_midpoint_fluid_sphcoord_mpt
