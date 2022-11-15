module interface_bh_boundary_psi_test2_mpt
  implicit none
  interface 
    subroutine bh_boundary_psi_test2_mpt(impt,char_bc,sou_surf,dsou_surf)
      integer          :: impt
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
      character(len=1), intent(in) :: char_bc
    end subroutine bh_boundary_psi_test2_mpt
  end interface
end module interface_bh_boundary_psi_test2_mpt
