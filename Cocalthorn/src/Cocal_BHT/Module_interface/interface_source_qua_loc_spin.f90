module interface_source_qua_loc_spin
  implicit none
  interface 
    subroutine source_qua_loc_spin(sous)
      real(8), pointer  :: sous(:,:)
    end subroutine source_qua_loc_spin
  end interface
end module interface_source_qua_loc_spin
