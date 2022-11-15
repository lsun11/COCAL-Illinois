module interface_poisson_solver_binary
  implicit none
  interface 
    subroutine poisson_solver_binary(sou,sou_surf,dsou_surf,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
    end subroutine poisson_solver_binary
  end interface
end module interface_poisson_solver_binary
