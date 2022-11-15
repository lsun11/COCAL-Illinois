module interface_poisson_solver_At
  implicit none
  interface 
    subroutine poisson_solver_At(sou,pot,char_sym)
      real(8),pointer  ::  sou(:,:,:), pot(:,:,:)
      character(len=6) ::  char_sym
    end subroutine poisson_solver_At
  end interface
end module interface_poisson_solver_At
