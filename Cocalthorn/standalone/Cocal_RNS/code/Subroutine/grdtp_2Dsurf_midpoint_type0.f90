subroutine grdtp_2Dsurf_midpoint_type0(fnc,deriv_theta,deriv_phi,itg,ipg)
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  implicit none
  real(8), pointer :: fnc(:,:)
  real(8), intent(out) :: deriv_theta, deriv_phi
  integer, intent(in)  :: itg, ipg
!
! --- Compute phi-derivative of a function.  
! --- The derivative is evaluated on grid points. 
!
  deriv_theta = 0.25d0 &
  &     *(fnc(itg  ,ipg) + fnc(itg  ,ipg-1) &
  &     + fnc(itg  ,ipg) + fnc(itg  ,ipg-1) &
  &     - fnc(itg-1,ipg) - fnc(itg-1,ipg-1) &
  &     - fnc(itg-1,ipg) - fnc(itg-1,ipg-1))*dthginv
!
  deriv_phi = 0.25d0 &
  &     *(fnc(itg  ,ipg) - fnc(itg  ,ipg-1) &
  &     + fnc(itg  ,ipg) - fnc(itg  ,ipg-1) &
  &     + fnc(itg-1,ipg) - fnc(itg-1,ipg-1) &
  &     + fnc(itg-1,ipg) - fnc(itg-1,ipg-1))*dphiginv
!
end subroutine grdtp_2Dsurf_midpoint_type0
