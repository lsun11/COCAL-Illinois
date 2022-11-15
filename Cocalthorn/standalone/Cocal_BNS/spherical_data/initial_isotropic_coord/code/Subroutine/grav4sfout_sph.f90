subroutine grav4sfout_sph(sousf,pot,iaho)
!
  use phys_constant, only : nnrg	!, pi
  use grid_parameter_1D, only : nrg
  implicit none
!
  real(8), intent(inout) :: pot(0:nnrg)
  real(8) :: sousf, pi4aho
  integer :: iaho
!
!  pi4inv = 1.0d0/4.0d0/pi
  pot(0:nrg) = pi4aho*pot(0:nrg)
!
end subroutine grav4sfout_sph
