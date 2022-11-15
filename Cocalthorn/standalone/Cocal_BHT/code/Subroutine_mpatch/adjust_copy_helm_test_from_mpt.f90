subroutine adjust_copy_helm_test_from_mpt(niq,msec_x_oold)
  use phys_constant,  only       : long, nmpt
  use grid_parameter, only       : rgin
  use def_binary_parameter, only : dis
  use def_bh_parameter, only     : ome_bh
  implicit none
  integer :: niq
  real(long) :: msec_x_oold(niq)
!
  call copy_def_binary_parameter_from_mpt(1)
  msec_x_oold(1) = dis
  call copy_grid_parameter_from_mpt(nmpt)
  call copy_def_binary_parameter_from_mpt(nmpt)
  call copy_def_bh_parameter_from_mpt(nmpt)
!
end subroutine adjust_copy_helm_test_from_mpt
