subroutine copy_weight_midpoint_fluid_to_mpt(impt)
  use phys_constant,  only : nnrf, nntf, nnpf
  use grid_parameter, only : nrf, ntf, npf
  use weight_midpoint_fluid
  use weight_midpoint_fluid_mpt
  use copy_array_static_1dto2d_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_arraystatic_1dto2d_mpt(impt, hwdrf, hwdrf_, 1, nnrf)
  call copy_arraystatic_1dto2d_mpt(impt, hwdtf, hwdtf_, 1, nntf)
  call copy_arraystatic_1dto2d_mpt(impt, hwdpf, hwdpf_, 1, nnpf)
  call copy_arraystatic_1dto2d_mpt(impt, tzwdrf, tzwdrf_, 0, nnrf)
  call copy_arraystatic_1dto2d_mpt(impt, siwdrf, siwdrf_, 0, nnrf)
  call copy_arraystatic_1dto2d_mpt(impt, siwdtf, siwdtf_, 0, nntf)
  call copy_arraystatic_1dto2d_mpt(impt, tzwdpf, tzwdpf_, 0, nnpf)
  call copy_arraystatic_1dto2d_mpt(impt, wdxf, wdxf_, 0, nnrf)
  call copy_array3dto4d_mpt(impt, hwrtpf, hwrtpf_, 1, nrf, 1, ntf, 1, npf)
  call copy_array3dto4d_mpt(impt, tzwrtpf, tzwrtpf_, 0, nrf, 1, ntf, 1, npf)
  call copy_array3dto4d_mpt(impt, siwrtpf, siwrtpf_, 0, nrf, 1, ntf, 1, npf)
  call copy_array3dto4d_mpt(impt,rtsiwrtpf,rtsiwrtpf_,0,nrf, 0, ntf, 0, npf)
!
end subroutine copy_weight_midpoint_fluid_to_mpt
