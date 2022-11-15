module interface_interpo_3D_initial_4th_surface
  implicit none
  interface 
    subroutine interpo_3D_initial_4th_surface(fnc)
      real(8), pointer :: fnc(:,:)
    end subroutine interpo_3D_initial_4th_surface
  end interface
end module interface_interpo_3D_initial_4th_surface
