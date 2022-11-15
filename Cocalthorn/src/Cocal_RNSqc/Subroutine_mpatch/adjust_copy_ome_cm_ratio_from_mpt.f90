subroutine adjust_copy_ome_cm_ratio_from_mpt(niq,msec_x_oold)
  use phys_constant,  only       : long, nmpt
  use grid_parameter, only       : rgin
  use def_binary_parameter, only : dis
  use def_bh_parameter, only     : ome_bh
  implicit none
  integer :: niq
  real(long) :: msec_x_oold(niq)
!
  call copy_def_bh_parameter_from_mpt(1)
  msec_x_oold(1) = ome_bh
  call copy_def_binary_parameter_from_mpt(1)
  msec_x_oold(2) = dis
  call copy_grid_parameter_from_mpt(2)
  msec_x_oold(3) = rgin
  call copy_grid_parameter_from_mpt(nmpt)
  call copy_def_binary_parameter_from_mpt(nmpt)
  call copy_def_bh_parameter_from_mpt(nmpt)
!
end subroutine adjust_copy_ome_cm_ratio_from_mpt
