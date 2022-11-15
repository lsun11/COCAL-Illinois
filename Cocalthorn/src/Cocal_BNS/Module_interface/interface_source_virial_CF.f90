module interface_source_virial_CF
  implicit none
  interface 
    subroutine source_virial_CF(sou_Tkin,sou_Pint,sou_Wgra)
      real(8), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
    end subroutine source_virial_CF
  end interface
end module interface_source_virial_CF
