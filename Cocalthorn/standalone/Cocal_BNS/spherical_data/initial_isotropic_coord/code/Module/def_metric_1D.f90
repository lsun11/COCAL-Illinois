module def_metric_1D
  use phys_constant, only : nnrg
  implicit none
!metpo
  real(8) :: psi(0:nnrg), alph(0:nnrg), alps(0:nnrg)	!alps_grav
!metpof
  real(8) :: psif(0:nnrg), alphf(0:nnrg)	!alps_fluid
end module def_metric_1D
