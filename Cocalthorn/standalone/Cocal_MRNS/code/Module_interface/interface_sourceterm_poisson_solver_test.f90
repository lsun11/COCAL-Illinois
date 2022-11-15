module interface_sourceterm_poisson_solver_test
  implicit none
  interface 
    subroutine sourceterm_poisson_solver_test(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_poisson_solver_test
  end interface
end module interface_sourceterm_poisson_solver_test
