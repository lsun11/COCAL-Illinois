module interface_grgrad_midpoint
  implicit none
  interface 
    subroutine grgrad_midpoint(fnc,dfdx,dfdy,dfdz)
      real(8), pointer :: fnc(:,:,:)
      real(8), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
    end subroutine grgrad_midpoint
  end interface
end module interface_grgrad_midpoint
