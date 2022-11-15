module interface_poisson_solver_binary_surf_int
  implicit none
  interface 
    subroutine poisson_solver_binary_surf_int(sou_surf,dsou_surf,pot)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:), pot(:,:,:)
    end subroutine poisson_solver_binary_surf_int
  end interface
end module interface_poisson_solver_binary_surf_int
