module interface_grgrad_gridpoint
  implicit none
  interface 
    subroutine grgrad_gridpoint(fnc,dfdx,dfdy,dfdz)
      real(8), pointer :: fnc(:,:,:)
      real(8), pointer :: dfdx(:,:,:)
      real(8), pointer :: dfdy(:,:,:)
      real(8), pointer :: dfdz(:,:,:)
    end subroutine grgrad_gridpoint
  end interface
end module interface_grgrad_gridpoint
