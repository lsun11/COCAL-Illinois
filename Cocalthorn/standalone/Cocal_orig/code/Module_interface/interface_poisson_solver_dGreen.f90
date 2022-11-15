module interface_poisson_solver_dGreen
  implicit none
  interface 
    subroutine poisson_solver_dGreen(ik, sou_bioa, pot)
      real(8),pointer  ::  sou_bioa(:,:,:,:), pot(:,:,:)
      integer :: ik
    end subroutine poisson_solver_dGreen
  end interface
end module interface_poisson_solver_dGreen
