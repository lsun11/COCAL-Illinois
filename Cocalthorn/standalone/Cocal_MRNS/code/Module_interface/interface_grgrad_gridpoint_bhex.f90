module interface_grgrad_gridpoint_bhex
  implicit none
  interface 
    subroutine grgrad_gridpoint_bhex(fnc,dfdx,dfdy,dfdz)
      real(8), pointer :: fnc(:,:,:)
      real(8), pointer :: dfdx(:,:,:)
      real(8), pointer :: dfdy(:,:,:)
      real(8), pointer :: dfdz(:,:,:)
    end subroutine grgrad_gridpoint_bhex
  end interface
end module interface_grgrad_gridpoint_bhex
