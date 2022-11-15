module interface_source_qua_loc_spin_peos_fluid
  implicit none
  interface 
    subroutine source_qua_loc_spin_peos_fluid(sous,irs)
      real(8), pointer  :: sous(:,:)
      integer :: irs
    end subroutine source_qua_loc_spin_peos_fluid
  end interface
end module interface_source_qua_loc_spin_peos_fluid
