module interface_source_qua_loc_spin_peos
  implicit none
  interface 
    subroutine source_qua_loc_spin_peos(sous,irs)
      real(8), pointer  :: sous(:,:)
      integer :: irs
    end subroutine source_qua_loc_spin_peos
  end interface
end module interface_source_qua_loc_spin_peos
