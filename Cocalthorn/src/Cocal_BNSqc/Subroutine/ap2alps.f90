subroutine ap2alps
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi, alph, alps, alps2
  implicit none
!
! --- Prepareing alps and alps2 from alpha and psi.  
!
  alps(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)
  alps2(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)**2
!
end subroutine ap2alps
