subroutine poisol(sousfo,dsousfo,soug,potg,iesy,irsy,imini,ieq)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrg
!  use coordinate_grav_r_1D, only : rg, hrg
  implicit none
!
  real(8), intent(inout) :: soug(0:nnrg),potg(0:nnrg)
  real(8), intent(inout) :: sousfo, dsousfo
  real(8) :: potg1(0:nnrg), poto1(0:nnrg),poto2(0:nnrg)
  integer :: iesy, irsy, imini, ieq
!
! --- Poisson solver : call volume and surface integrals.
!
!     sources for volume integral are on half integer points.
!     sources for surface integral are on half integer points.
!
! --- Call integral terms of Green's formula.
!
!  call halfsou(soug)
  call interpo_fl2grmidpoint_1D(soug)
!
  call gravmid_sph(soug,potg1)
  call grav4sfout_sph(dsousfo,poto1,3)
  call grav4sfout_sph( sousfo,poto2,4)
!
  potg(0:nrg) = potg1(0:nrg) + 1.0d0
!testtest     &          + poto1(0:nrg) + poto2(0:nrg)
!
end subroutine poisol
