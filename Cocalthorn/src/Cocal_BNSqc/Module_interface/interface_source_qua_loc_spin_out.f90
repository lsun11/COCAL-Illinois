module interface_source_qua_loc_spin_out
  implicit none
  interface 
    subroutine source_qua_loc_spin_out(sous,irs)
      real(8), pointer  :: sous(:,:)
      integer :: irs
    end subroutine source_qua_loc_spin_out
  end interface
end module interface_source_qua_loc_spin_out
