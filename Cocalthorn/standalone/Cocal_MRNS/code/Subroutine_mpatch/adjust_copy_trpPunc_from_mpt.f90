subroutine adjust_copy_trpPunc_from_mpt(niq,msec_x_oold)
  use phys_constant,  only       : long, nmpt
  use def_bh_parameter, only     : mom_pBH, mass_pBH
  use def_binary_parameter, only : dis
  implicit none
  integer :: niq
  real(long) :: msec_x_oold(niq)
!
  call copy_def_bh_parameter_from_mpt(1)
  msec_x_oold(1) = mom_pBH(2)
  call copy_def_binary_parameter_from_mpt(1)
  msec_x_oold(2) = dis
!  call copy_def_bh_parameter_from_mpt(2)
!  msec_x_oold(2) = mom_pBH(2)
  call copy_def_bh_parameter_from_mpt(2)
  msec_x_oold(3) = mass_pBH
  call copy_def_bh_parameter_from_mpt(nmpt)
  call copy_def_binary_parameter_from_mpt(nmpt)
  call copy_def_quantities_bh_from_mpt(nmpt)
!
end subroutine adjust_copy_trpPunc_from_mpt
