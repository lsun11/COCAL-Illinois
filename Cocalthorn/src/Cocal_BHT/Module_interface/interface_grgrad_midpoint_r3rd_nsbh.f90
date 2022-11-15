module interface_grgrad_midpoint_r3rd_nsbh
  implicit none
  interface 
    subroutine grgrad_midpoint_r3rd_nsbh(fnc,dfdx,dfdy,dfdz,cobj)
      real(8), pointer :: fnc(:,:,:)
      real(8), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine grgrad_midpoint_r3rd_nsbh
  end interface
end module interface_grgrad_midpoint_r3rd_nsbh
