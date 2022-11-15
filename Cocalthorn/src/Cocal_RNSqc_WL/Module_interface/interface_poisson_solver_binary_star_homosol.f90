module interface_poisson_solver_binary_star_homosol
  implicit none
  interface 
    subroutine poisson_solver_binary_star_homosol(char_bc,sou,&
                                            &     sou_exsurf,dsou_exsurf, &
                                            &    sou_outsurf,dsou_outsurf,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
      real(8), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
      character(len=2) :: char_bc
    end subroutine poisson_solver_binary_star_homosol
  end interface
end module interface_poisson_solver_binary_star_homosol
