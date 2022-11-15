subroutine adjust_calc_fncval_helm_test_mpt(niq,msec_f_oold)
  use phys_constant, only  : long, nmpt
  use def_quantities, only : admmom_asymp
  implicit none
  integer :: niq
  real(long) :: msec_f_oold(niq)
!
  call copy_def_quantities_from_mpt(nmpt)
!xpsi  msec_f_oold(1) = admmom_asymp(1)
  msec_f_oold(1) = admmom_asymp(2)
!
end subroutine adjust_calc_fncval_helm_test_mpt
