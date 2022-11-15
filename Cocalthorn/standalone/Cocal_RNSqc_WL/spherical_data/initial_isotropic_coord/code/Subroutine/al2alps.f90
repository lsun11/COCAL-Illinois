subroutine al2alps
!
  use def_metric_1D, only : alps, alph, psi
  use grid_parameter_1D, only : nrg
  implicit none
!
! --  set alps
!
  alps(0:nrg) = alph(0:nrg)*psi(0:nrg)
!
end subroutine al2alps
