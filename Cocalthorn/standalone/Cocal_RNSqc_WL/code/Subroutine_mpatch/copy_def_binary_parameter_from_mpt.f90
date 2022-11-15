subroutine copy_def_binary_parameter_from_mpt(impt)
  use def_binary_parameter
  use def_binary_parameter_mpt
  use copy_array_static_1dto0d_mpt
  implicit none
  integer :: impt
!
  call copy_arraystatic_1dto0d_mpt(impt, sepa_, sepa)
  call copy_arraystatic_1dto0d_mpt(impt, dis_, dis)
  call copy_arraystatic_1dto0d_mpt(impt, mass_ratio_, mass_ratio)
!
end subroutine copy_def_binary_parameter_from_mpt
