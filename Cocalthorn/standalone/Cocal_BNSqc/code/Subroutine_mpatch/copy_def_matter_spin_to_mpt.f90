subroutine copy_def_matter_spin_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg, nlg, nrf, ntf, npf
  use def_metric
  use def_metric_mpt
  use def_matter
  use def_matter_mpt
  use def_velocity_potential
  use def_velocity_potential_mpt
  use def_velocity_rot
  use def_velocity_rot_mpt
  use copy_array_4dto5d_mpt
  use copy_array_3dto4d_mpt
  use copy_array_2dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt, vep  , vep_  , 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, vepxf, vepxf_, 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, vepyf, vepyf_, 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, vepzf, vepzf_, 0, nrf, 0, ntf, 0, npf)
!
  call copy_array3dto4d_mpt(impt, vepxg, vepxg_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, vepyg, vepyg_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, vepzg, vepzg_, 0, nrg, 0, ntg, 0, npg)

  call copy_array3dto4d_mpt(impt, wxspf, wxspf_, 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, wyspf, wyspf_, 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, wzspf, wzspf_, 0, nrf, 0, ntf, 0, npf)
!
  call copy_array3dto4d_mpt(impt, wxspg, wxspg_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, wyspg, wyspg_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, wzspg, wzspg_, 0, nrg, 0, ntg, 0, npg)

end subroutine copy_def_matter_spin_to_mpt
