subroutine grdphi_midpoint_type0(fnc,deriv,irg,itg,ipg)
  use coordinate_grav_phi, only : dphiginv
  implicit none
  real(8), pointer :: fnc(:,:,:)
  real(8), intent(out) :: deriv
  integer, intent(in) :: irg, itg, ipg
!
! --- Compute phi-derivative of a function.  
! --- The derivative is evaluated on grid points. 
!
  deriv = 0.25d0 &
  &     *(fnc(irg  ,itg  ,ipg) - fnc(irg  ,itg  ,ipg-1) &
  &     + fnc(irg-1,itg  ,ipg) - fnc(irg-1,itg  ,ipg-1) &
  &     + fnc(irg  ,itg-1,ipg) - fnc(irg  ,itg-1,ipg-1) &
  &     + fnc(irg-1,itg-1,ipg) - fnc(irg-1,itg-1,ipg-1))*dphiginv
!
end subroutine grdphi_midpoint_type0
