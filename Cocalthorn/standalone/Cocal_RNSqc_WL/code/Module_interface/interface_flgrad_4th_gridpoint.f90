module interface_flgrad_4th_gridpoint
  implicit none
  interface
    subroutine flgrad_4th_gridpoint(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irf, itf, ipf
    end subroutine flgrad_4th_gridpoint
  end interface
end module interface_flgrad_4th_gridpoint
