module interface_outer_boundary_d_poty
  implicit none
  interface 
    subroutine outer_boundary_d_poty(sou_surf,Bfun,dBfundy)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:), Bfun(:,:,:), dBfundy(:,:,:)
    end subroutine outer_boundary_d_poty
  end interface
end module interface_outer_boundary_d_poty
