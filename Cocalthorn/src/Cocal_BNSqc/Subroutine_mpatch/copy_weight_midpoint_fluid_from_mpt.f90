subroutine copy_weight_midpoint_fluid_from_mpt(impt)
  use phys_constant,  only : nnrf, nntf, nnpf
  use grid_parameter, only : nrf, ntf, npf
  use weight_midpoint_fluid
  use weight_midpoint_fluid_mpt
  use copy_array_static_2dto1d_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_arraystatic_2dto1d_mpt(impt, hwdrf_, hwdrf, 1, nnrf)
  call copy_arraystatic_2dto1d_mpt(impt, hwdtf_, hwdtf, 1, nntf)
  call copy_arraystatic_2dto1d_mpt(impt, hwdpf_, hwdpf, 1, nnpf)
  call copy_arraystatic_2dto1d_mpt(impt, tzwdrf_, tzwdrf, 0, nnrf)
  call copy_arraystatic_2dto1d_mpt(impt, siwdrf_, siwdrf, 0, nnrf)
  call copy_arraystatic_2dto1d_mpt(impt, siwdtf_, siwdtf, 0, nntf)
  call copy_arraystatic_2dto1d_mpt(impt, tzwdpf_, tzwdpf, 0, nnpf)
  call copy_arraystatic_2dto1d_mpt(impt, wdxf_, wdxf, 0, nnrf)
  call copy_array4dto3d_mpt(impt, hwrtpf_, hwrtpf, 1, nrf, 1, ntf, 1, npf)
  call copy_array4dto3d_mpt(impt, tzwrtpf_, tzwrtpf, 0, nrf, 1, ntf, 1, npf)
  call copy_array4dto3d_mpt(impt, siwrtpf_, siwrtpf, 0, nrf, 1, ntf, 1, npf)
  call copy_array4dto3d_mpt(impt, rtsiwrtpf_,rtsiwrtpf, 0, nrf, 0, ntf, 0, npf)
!
end subroutine copy_weight_midpoint_fluid_from_mpt
