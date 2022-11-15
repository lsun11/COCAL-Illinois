module interface_source_virial_WL
  implicit none
  interface 
    subroutine source_virial_WL(sou_Tkin,sou_Pint,sou_Wgra)
      real(8), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
    end subroutine source_virial_WL
  end interface
end module interface_source_virial_WL
