module interface_outer_boundary_d_BBH_CF_pBH
  implicit none
  interface 
    subroutine outer_boundary_d_BBH_CF_pBH(sou_surf,char_mp)
      real(8), pointer :: sou_surf(:,:)
      character(len=4), intent(in) :: char_mp
    end subroutine outer_boundary_d_BBH_CF_pBH
  end interface
end module interface_outer_boundary_d_BBH_CF_pBH
