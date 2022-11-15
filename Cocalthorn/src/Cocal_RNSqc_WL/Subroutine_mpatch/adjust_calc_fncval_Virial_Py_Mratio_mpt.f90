subroutine adjust_calc_fncval_Virial_Py_Mratio_mpt(niq,msec_f_oold)
  use phys_constant, only  : long, nmpt
  use def_quantities, only : admmass_asymp, komarmass_asymp, admmom_asymp
  use def_quantities_bh,    only : AHmass
  use def_binary_parameter, only : mass_ratio
  implicit none
  integer :: niq
  real(long) :: msec_f_oold(niq), Mirr1, Mirr2
!
  call copy_def_quantities_from_mpt(nmpt)
  msec_f_oold(1) = (admmass_asymp - komarmass_asymp)/admmass_asymp             
  msec_f_oold(2) = admmom_asymp(2)
  call copy_def_quantities_bh_from_mpt(1)
  Mirr1 = AHmass
  call copy_def_quantities_bh_from_mpt(2)
  Mirr2 = AHmass
  msec_f_oold(3) = Mirr2/(Mirr1*mass_ratio) - 1.0d0
!
end subroutine adjust_calc_fncval_Virial_Py_Mratio_mpt
