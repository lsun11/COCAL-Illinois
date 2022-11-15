module interface_bh_boundary_n_Bfun
  implicit none
  interface 
    subroutine bh_boundary_n_Bfun(dsou_surf,potx,poty,potz)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: dsou_surf(:,:), potx(:,:,:), poty(:,:,:), potz(:,:,:)
    end subroutine bh_boundary_n_Bfun
  end interface
end module interface_bh_boundary_n_Bfun
