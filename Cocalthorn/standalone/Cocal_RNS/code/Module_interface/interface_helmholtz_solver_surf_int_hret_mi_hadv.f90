module interface_helmholtz_solver_surf_int_hret_mi_hadv
  implicit none
  interface 
    subroutine helmholtz_solver_surf_int_hret_mi_hadv(sou_surf,dsou_surf,pot)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:), pot(:,:,:)
    end subroutine helmholtz_solver_surf_int_hret_mi_hadv
  end interface
end module interface_helmholtz_solver_surf_int_hret_mi_hadv
