subroutine compute_psi2wme
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,     only : psi, alph
  use def_metric_pBH, only : wme, log_wme, log_N, index_wme
  use interface_compute_fnc_inversion_index
  implicit none
  real(long) :: index, fac2
  integer :: irg, itg, ipg
!
  index = - dble(index_wme)
  call compute_fnc_inversion_index(psi,wme,index)
  wme(0,0:ntg,0:npg) = 0.0d0
  log_wme(1:nrg,0:ntg,0:npg) = dlog(wme(1:nrg,0:ntg,0:npg))
  log_wme(0,    0:ntg,0:npg) = -40.0d0
  fac2 = 1.0d0/dsqrt(2.0d0)
  log_N(1:nrg,0:ntg,0:npg) = fac2*dlog(alph(1:nrg,0:ntg,0:npg))
  log_N(0,    0:ntg,0:npg) = -40.0d0
!
end subroutine compute_psi2wme
