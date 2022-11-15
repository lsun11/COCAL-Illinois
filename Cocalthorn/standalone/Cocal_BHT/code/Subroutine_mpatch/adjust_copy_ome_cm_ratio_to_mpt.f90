subroutine adjust_copy_ome_cm_ratio_to_mpt(niq,msec_x)
  use phys_constant, only  : long, nmpt
  use grid_parameter, only : rgin
  use def_binary_parameter, only : dis, sepa
  use def_bh_parameter, only : ome_bh
  implicit none
  integer :: niq, ii
  real(long) :: msec_x(niq)
!
  do ii = 1, nmpt
    call copy_grid_parameter_from_mpt(ii)
    call copy_def_binary_parameter_from_mpt(ii)
    call copy_def_bh_parameter_from_mpt(ii)
    if (ii.eq.1) then 
      ome_bh = msec_x(1)
      dis    = msec_x(2)
    end if
    if (ii.eq.2) then 
      ome_bh = msec_x(1)
      dis    = sepa - msec_x(2)
      rgin   = msec_x(3)
    end if
    if (ii.eq.nmpt) then 
      ome_bh = msec_x(1)
    end if
    call copy_grid_parameter_to_mpt(ii)
    call copy_def_binary_parameter_to_mpt(ii)
    call copy_def_bh_parameter_to_mpt(ii)
  end do
!
end subroutine adjust_copy_ome_cm_ratio_to_mpt
