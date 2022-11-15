module interface_bh_boundary_d_potx
  implicit none
  interface 
    subroutine bh_boundary_d_potx(sou_surf,Bfun,dBfundx)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:), Bfun(:,:,:), dBfundx(:,:,:)
    end subroutine bh_boundary_d_potx
  end interface
end module interface_bh_boundary_d_potx
