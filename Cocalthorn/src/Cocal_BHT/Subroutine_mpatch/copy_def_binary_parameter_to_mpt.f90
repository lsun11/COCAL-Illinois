subroutine copy_def_binary_parameter_to_mpt(impt)
  use def_binary_parameter
  use def_binary_parameter_mpt
  use copy_array_static_0dto1d_mpt
  implicit none
  integer :: impt
!
  call copy_arraystatic_0dto1d_mpt(impt, sepa, sepa_)
  call copy_arraystatic_0dto1d_mpt(impt, dis, dis_)
  call copy_arraystatic_0dto1d_mpt(impt, mass_ratio, mass_ratio_)
!
end subroutine copy_def_binary_parameter_to_mpt
