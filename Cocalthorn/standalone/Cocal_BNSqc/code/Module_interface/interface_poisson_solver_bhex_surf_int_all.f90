module interface_poisson_solver_bhex_surf_int_all
  implicit none
  interface 
    subroutine poisson_solver_bhex_surf_int_all(sou_iosurf,pot)
      real(8), pointer :: sou_iosurf(:,:,:), pot(:,:,:)
    end subroutine poisson_solver_bhex_surf_int_all
  end interface
end module interface_poisson_solver_bhex_surf_int_all
