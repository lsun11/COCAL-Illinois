module interface_poisson_solver_binary_bhex_homosol
  implicit none
  interface 
    subroutine poisson_solver_binary_bhex_homosol(char_bc,sou,&
                                            &     sou_exsurf,dsou_exsurf, &
                                            &     sou_insurf,dsou_insurf, &
                                            &    sou_outsurf,dsou_outsurf,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
      real(8), pointer :: sou_insurf(:,:), dsou_insurf(:,:)
      real(8), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
      character(len=2) :: char_bc
    end subroutine poisson_solver_binary_bhex_homosol
  end interface
end module interface_poisson_solver_binary_bhex_homosol
