module interface_poisson_solver_homogeneous_sol_spin
  implicit none
  interface 
    subroutine poisson_solver_homogeneous_sol_spin(surp, vpot_b)
      real(8),pointer  ::  surp(:,:), vpot_b(:,:,:)
    end subroutine  poisson_solver_homogeneous_sol_spin
  end interface
end module interface_poisson_solver_homogeneous_sol_spin
