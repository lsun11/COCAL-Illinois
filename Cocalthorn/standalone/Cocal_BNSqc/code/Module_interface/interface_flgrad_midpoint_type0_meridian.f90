module interface_flgrad_midpoint_type0_meridian
  implicit none
  interface 
    subroutine flgrad_midpoint_type0_meridian(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
      real(8), pointer :: fnc(:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irf, itf, ipf
    end subroutine flgrad_midpoint_type0_meridian
  end interface
end module interface_flgrad_midpoint_type0_meridian
