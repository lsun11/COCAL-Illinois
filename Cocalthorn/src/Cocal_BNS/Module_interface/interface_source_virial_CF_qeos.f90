module interface_source_virial_CF_qeos
  implicit none
  interface 
    subroutine source_virial_CF_qeos(sou_Tkin,sou_Pint,sou_Wgra)
      real(8), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
    end subroutine source_virial_CF_qeos
  end interface
end module interface_source_virial_CF_qeos
