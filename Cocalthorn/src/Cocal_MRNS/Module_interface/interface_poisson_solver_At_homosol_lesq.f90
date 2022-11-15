module interface_poisson_solver_At_homosol_lesq
  implicit none
  interface 
    subroutine poisson_solver_At_homosol_lesq(pot_bc,pot_nb,pot,char_sym)
      real(8),pointer  ::  pot_bc(:,:), pot_nb(:,:), pot(:,:,:)
      character(len=6) ::  char_sym
    end subroutine poisson_solver_At_homosol_lesq
  end interface
end module interface_poisson_solver_At_homosol_lesq
