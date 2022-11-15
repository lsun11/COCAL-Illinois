module interface_poisson_solver_1bh_hij
  implicit none
  interface 
    subroutine poisson_solver_1bh_hij(ik, char_bc, sou, sou_bioa,  &
                                &     sou_insurf, dsou_insurf, &
                                &     sou_outsurf, dsou_outsurf, dsou_inbaij, pot)
      real(8), pointer :: sou_bioa(:,:,:,:)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_insurf(:,:), dsou_insurf(:,:), dsou_inbaij(:,:)
      real(8), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
      character(len=2) :: char_bc
      integer :: ik
    end subroutine poisson_solver_1bh_hij
  end interface
end module interface_poisson_solver_1bh_hij
