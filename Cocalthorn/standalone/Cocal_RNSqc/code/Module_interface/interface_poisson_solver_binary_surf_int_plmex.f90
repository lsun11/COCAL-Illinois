module interface_poisson_solver_binary_surf_int_plmex
  implicit none
  interface 
    subroutine poisson_solver_binary_surf_int_plmex(impt,sou_surf,dsou_surf,pot)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:), pot(:,:,:)
      integer :: impt
    end subroutine poisson_solver_binary_surf_int_plmex
  end interface
end module interface_poisson_solver_binary_surf_int_plmex
