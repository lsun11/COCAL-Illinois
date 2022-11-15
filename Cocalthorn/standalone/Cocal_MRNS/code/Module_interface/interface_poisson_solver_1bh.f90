module interface_poisson_solver_1bh
  implicit none
  interface 
    subroutine poisson_solver_1bh(char_bc, sou,&
                            &     sou_insurf, dsou_insurf, &
                            &     sou_outsurf, dsou_outsurf, pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_insurf(:,:), dsou_insurf(:,:)
      real(8), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
      character(len=2) :: char_bc
    end subroutine poisson_solver_1bh
  end interface
end module interface_poisson_solver_1bh
