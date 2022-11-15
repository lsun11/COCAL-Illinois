module interface_poisson_solver_dGreen_v1
  implicit none
  interface 
    subroutine poisson_solver_dGreen_v1(ik, sou_bioa, pot)
      real(8),pointer  ::  sou_bioa(:,:,:,:), pot(:,:,:)
      integer :: ik
    end subroutine poisson_solver_dGreen_v1
  end interface
end module interface_poisson_solver_dGreen_v1
