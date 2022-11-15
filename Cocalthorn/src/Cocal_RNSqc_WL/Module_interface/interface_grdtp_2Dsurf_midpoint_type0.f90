module interface_grdtp_2Dsurf_midpoint_type0
  implicit none
  interface 
    subroutine grdtp_2Dsurf_midpoint_type0(fnc,deriv_theta,deriv_phi,itg,ipg)
      real(8), pointer :: fnc(:,:)
      real(8), intent(out) :: deriv_theta, deriv_phi
      integer, intent(in)  :: itg, ipg
    end subroutine grdtp_2Dsurf_midpoint_type0
  end interface
end module interface_grdtp_2Dsurf_midpoint_type0
