module interface_source_virial_gravity_WL
  implicit none
  interface 
    subroutine source_virial_gravity_WL(sou_Wgra)
      real(8), pointer :: sou_Wgra(:,:,:)
    end subroutine source_virial_gravity_WL
  end interface
end module interface_source_virial_gravity_WL
