subroutine copy_def_matter_irrot_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg, nlg, nrf, ntf, npf
  use def_metric
  use def_metric_mpt
  use def_matter
  use def_matter_mpt
  use def_velocity_potential
  use def_velocity_potential_mpt
  use copy_array_5dto4d_mpt
  use copy_array_4dto3d_mpt
  use copy_array_3dto2d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, vep_  , vep  , 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, vepxf_, vepxf, 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, vepyf_, vepyf, 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, vepzf_, vepzf, 0, nrf, 0, ntf, 0, npf)
!
  call copy_array4dto3d_mpt(impt, vepxg_, vepxg, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, vepyg_, vepyg, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, vepzg_, vepzg, 0, nrg, 0, ntg, 0, npg)

end subroutine copy_def_matter_irrot_from_mpt
