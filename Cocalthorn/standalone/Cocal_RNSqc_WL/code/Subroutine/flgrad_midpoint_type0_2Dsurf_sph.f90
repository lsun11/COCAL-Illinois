subroutine flgrad_midpoint_type0_2Dsurf_sph(fnc,dfdth,dfdphi,itf,ipf)
  use phys_constant, only : long
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  implicit none
  real(long), pointer :: fnc(:,:)
  real(long), intent(out) :: dfdth, dfdphi
  integer, intent(in) :: itf, ipf
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
  dfdth  = 0.5d0 &
  &       *(fnc(itf,ipf  ) - fnc(itf-1,ipf  ) &
  &       + fnc(itf,ipf-1) - fnc(itf-1,ipf-1))*dthginv
  dfdphi = 0.5d0 &
  &       *(fnc(itf  ,ipf) - fnc(itf  ,ipf-1) &
  &       + fnc(itf-1,ipf) - fnc(itf-1,ipf-1))*dphiginv
!
end subroutine flgrad_midpoint_type0_2Dsurf_sph
