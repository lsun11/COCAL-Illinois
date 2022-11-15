subroutine calc_physq_center_peos
  use def_matter, only : emd
  use def_quantities, only : rho_c, pre_c, epsi_c, q_c 
  implicit none
  real(8) :: hh, emdmax, xsol
  integer    :: ir, iremax
!
  call search_emdmax_xaxis_grid(iremax)
  if (iremax.eq.0) then 
    q_c = emd(0,0,0)
  else
    call search_emdmax_xaxis(iremax,emdmax,xsol)
    q_c = emdmax
  end if
  call peos_q2hprho(q_c, hh, pre_c, rho_c, epsi_c)
!
end subroutine calc_physq_center_peos
