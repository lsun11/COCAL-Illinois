module interface_poisson_solver_AHfinder
  implicit none
  interface 
    subroutine poisson_solver_AHfinder(sou,pot)
      real(8),pointer  ::  sou(:,:), pot(:,:)
    end subroutine poisson_solver_AHfinder
  end interface
end module interface_poisson_solver_AHfinder
