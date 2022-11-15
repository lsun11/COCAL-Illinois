module interface_source_qua_loc_spin_grav
  implicit none
  interface 
    subroutine source_qua_loc_spin_grav(sous,irs)
      real(8), pointer  :: sous(:,:)
      integer :: irs
    end subroutine source_qua_loc_spin_grav
  end interface
end module interface_source_qua_loc_spin_grav
