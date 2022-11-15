module interface_grgrad_midpoint_r3rd_type0_bh
  implicit none
  interface 
    subroutine grgrad_midpoint_r3rd_type0_bh(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_midpoint_r3rd_type0_bh
  end interface
end module interface_grgrad_midpoint_r3rd_type0_bh
