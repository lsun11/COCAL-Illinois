module def_qeos_parameter
  use phys_constant         !nnpeos
  implicit none
  real(8) :: abc(0:nnqeos), abi(0:nnqeos), rhoi(0:nnqeos), &
  &          qi(0:nnqeos), hi(0:nnqeos)
  real(8) :: abccgs(0:nnqeos)
  real(8) :: abcene(0:nnqeos), abiene(0:nnqeos), abch(0:nnqeos), abih(0:nnqeos)
  real(8) :: abchdot(0:nnqeos), abihdot(0:nnqeos)
  real(8) :: rhoini_cgs, rhoini_gcm1, emdini_gcm1  !used in TOV solver
  real(8) :: eneini_cgs, eneini_gcm1 !used in MIT bag model
  real(8) :: rhosurf_cgs, rhosurf_gcm1, hsurf_gcm1 ! used for RNS_qeos
  real(8) :: eneconst_cgs, eneconst_gcm1
  real(8) :: enesurf_cgs, enesurf_gcm1 ! used in MIT bag model
  real(8) :: aq                        ! used in MIT bag model
  integer :: nphase
end module def_qeos_parameter
