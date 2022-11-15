module interface_bh_boundary_d_potz
  implicit none
  interface 
    subroutine bh_boundary_d_potz(sou_surf,Bfun,dBfundz)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:), Bfun(:,:,:), dBfundz(:,:,:)
    end subroutine bh_boundary_d_potz
  end interface
end module interface_bh_boundary_d_potz
