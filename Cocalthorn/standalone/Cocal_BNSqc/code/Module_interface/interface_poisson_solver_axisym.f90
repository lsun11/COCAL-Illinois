module interface_poisson_solver_axisym
  implicit none
  interface 
    subroutine poisson_solver_axisym(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver_axisym
  end interface
end module interface_poisson_solver_axisym
