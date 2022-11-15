module interface_flgrad_midpoint_surface_type0
  implicit none
  interface 
    subroutine flgrad_midpoint_surface_type0(fnc,dfdx,dfdy,dfdz,itf,ipf)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: itf, ipf
    end subroutine flgrad_midpoint_surface_type0
  end interface
end module interface_flgrad_midpoint_surface_type0
