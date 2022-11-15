!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_binary_excision_mpt
  use phys_constant, only : long
  implicit none
! weight binary excision
  real(long), pointer :: hwrtpg_ex_(:,:,:,:)
  real(long), pointer :: hwtpg_ex_(:,:,:)
end module weight_midpoint_binary_excision_mpt
