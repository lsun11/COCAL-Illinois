module interface_source_virial_gravity_CF
  implicit none
  interface 
    subroutine source_virial_gravity_CF(sou_Wgra)
      real(8), pointer :: sou_Wgra(:,:,:)
    end subroutine source_virial_gravity_CF
  end interface
end module interface_source_virial_gravity_CF
