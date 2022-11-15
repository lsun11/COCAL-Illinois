module interface_grgrad_midpoint_r2nd
  implicit none
  interface 
    subroutine grgrad_midpoint_r2nd(fnc,dfdx,dfdy,dfdz)
      real(8), pointer :: fnc(:,:,:)
      real(8), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
    end subroutine grgrad_midpoint_r2nd
  end interface
end module interface_grgrad_midpoint_r2nd
