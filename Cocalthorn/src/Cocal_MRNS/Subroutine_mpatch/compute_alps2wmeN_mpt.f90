subroutine compute_alps2wmeN_mpt(impt)
  use phys_constant,  only : long, nmpt
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,     only : psi, alph
  use def_metric_pBH, only : wme, log_wme, log_N, index_wme
  use interface_compute_fnc_inversion_index
  implicit none
  real(long) :: index, fac2
  integer :: irg, itg, ipg, impt
!
  index = - dble(index_wme)
  fac2 = 1.0d0/dsqrt(2.0d0)
  call compute_fnc_inversion_index(psi,wme,index)
  log_wme(1:nrg,0:ntg,0:npg) = dlog(wme(1:nrg,0:ntg,0:npg))
  log_N(  1:nrg,0:ntg,0:npg) = fac2*dlog(alph(1:nrg,0:ntg,0:npg))
  if (impt.ne.nmpt) then
        wme(0,0:ntg,0:npg) = 0.0d0
    log_wme(0,0:ntg,0:npg) = -40.0d0
    log_N(  0,0:ntg,0:npg) = -40.0d0
  else 
    log_wme(0,0:ntg,0:npg) = dlog(wme(0,0:ntg,0:npg))
    log_N(  0,0:ntg,0:npg) = fac2*dlog(alph(0,0:ntg,0:npg))
  end if
!
end subroutine compute_alps2wmeN_mpt
