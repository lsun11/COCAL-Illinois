subroutine alps2al
!
  use def_metric_1D, only : alph, alps, psi
  use grid_parameter_1D, only : nrg
  implicit none
!
! --  set alpha
!
    alph(0:nrg) = alps(0:nrg)/psi(0:nrg)
!   alph      alps(irg) = alph(irg)*psi(irg)
!
end subroutine alps2al
