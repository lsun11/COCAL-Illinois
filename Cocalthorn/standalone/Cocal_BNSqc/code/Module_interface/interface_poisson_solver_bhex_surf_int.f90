module interface_poisson_solver_bhex_surf_int
  implicit none
  interface 
    subroutine poisson_solver_bhex_surf_int(char_io,sou_iosurf,pot)
      real(8), pointer :: sou_iosurf(:,:,:), pot(:,:,:)
      character(len=2), intent(in) :: char_io
    end subroutine poisson_solver_bhex_surf_int
  end interface
end module interface_poisson_solver_bhex_surf_int
