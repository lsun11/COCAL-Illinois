subroutine copy_def_matter_to_mpt(impt)
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_matter
  use def_matter_mpt
  use def_velocity_rot
  use def_velocity_rot_mpt
  use copy_array_3dto4d_mpt
  use copy_array_2dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array2dto3d_mpt(impt, rs,        rs_,        0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, emd,       emd_,       0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, utf,       utf_,       0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, omef,      omef_,      0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, jomef,     jomef_,     0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, jomef_int, jomef_int_, 0, nrf, 0, ntf, 0, npf)

  call copy_array3dto4d_mpt(impt, emdg,      emdg_,      0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, utg,       utg_,       0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, omeg,      omeg_,      0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, jomeg,     jomeg_,     0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, jomeg_int, jomeg_int_, 0, nrg, 0, ntg, 0, npg)
!
  call copy_array3dto4d_mpt(impt, vep,       vep_,       0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, wxspf,   wxspf_,       0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, wyspf,   wyspf_,       0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, wzspf,   wzspf_,       0, nrf, 0, ntf, 0, npf)
!
end subroutine copy_def_matter_to_mpt
