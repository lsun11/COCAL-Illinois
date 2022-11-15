subroutine adjust_copy_trpPunc_to_mpt(niq,msec_x)
  use phys_constant, only  : long, nmpt
  use def_bh_parameter,  only : mom_pBH, mass_pBH
  use def_binary_parameter, only : dis, sepa
  implicit none
  integer :: niq, ii
  real(long) :: msec_x(niq)
!
  do ii = 1, 2
    call copy_def_bh_parameter_from_mpt(ii)
    call copy_def_binary_parameter_from_mpt(ii)
    if (ii.eq.1) then 
      mom_pBH(2) = msec_x(1)
      dis    = msec_x(2)
    end if
    if (ii.eq.2) then 
      mom_pBH(2) = msec_x(1)
      dis        = sepa - msec_x(2)
      mass_pBH   = msec_x(3)
    end if
    call copy_def_bh_parameter_to_mpt(ii)
    call copy_def_binary_parameter_to_mpt(ii)
  end do
!
end subroutine adjust_copy_trpPunc_to_mpt
