module interface_helmholtz_solver_outer_surf_int
  implicit none
  interface 
    subroutine helmholtz_solver_outer_surf_int(sou_surf,dsou_surf,pot)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:), pot(:,:,:)
    end subroutine helmholtz_solver_outer_surf_int
  end interface
end module interface_helmholtz_solver_outer_surf_int
