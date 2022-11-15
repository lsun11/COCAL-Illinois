module interface_flgrad_2nd_gridpoint_export
  implicit none
  interface 
    subroutine flgrad_2nd_gridpoint_export(fnc,dfdx,dfdy,dfdz,irf,itf,ipf,rs)
      real(8), pointer :: fnc(:,:,:), rs(:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irf, itf, ipf
    end subroutine flgrad_2nd_gridpoint_export
  end interface
end module interface_flgrad_2nd_gridpoint_export
