module interface_flgrad_2nd_gridpoint
  implicit none
  interface 
    subroutine flgrad_2nd_gridpoint(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irf, itf, ipf
    end subroutine flgrad_2nd_gridpoint
  end interface
end module interface_flgrad_2nd_gridpoint
