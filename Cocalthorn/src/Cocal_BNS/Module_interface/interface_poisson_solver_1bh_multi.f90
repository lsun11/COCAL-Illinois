module interface_poisson_solver_1bh_multi
  implicit none
  interface 
    subroutine poisson_solver_1bh_multi(mchar_bc, msou,&
                                  &     msou_insurf, mdsou_insurf, &
                                  &     msou_outsurf, mdsou_outsurf, mpot)
      real(8), pointer :: mpot(:,:,:,:), msou(:,:,:,:)
      real(8), pointer :: msou_insurf(:,:,:), mdsou_insurf(:,:,:)
      real(8), pointer :: msou_outsurf(:,:,:), mdsou_outsurf(:,:,:)
      character(len=2), pointer :: mchar_bc(:)
    end subroutine poisson_solver_1bh_multi
  end interface
end module interface_poisson_solver_1bh_multi
