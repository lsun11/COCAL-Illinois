!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_fluid_sphcoord
  use phys_constant, only : long
  implicit none
! weight binary excision
  real(long), pointer :: hwrtpg_fc(:,:,:)
  real(long), pointer :: hwtpg_fc(:,:)
!
contains
!
subroutine allocate_weight_midpoint_fluid_sphcoord
  use grid_parameter, only : nrf, ntg, npg 
  use make_array_2d
  use make_array_3d
  implicit none
  call alloc_array3d(hwrtpg_fc,1,nrf,1,ntg,1,npg)
  call alloc_array2d(hwtpg_fc,1,ntg,1,npg)
end subroutine allocate_weight_midpoint_fluid_sphcoord
!
end module weight_midpoint_fluid_sphcoord
