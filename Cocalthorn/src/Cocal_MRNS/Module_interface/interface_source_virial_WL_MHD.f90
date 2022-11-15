module interface_source_virial_WL_MHD
  implicit none
  interface 
    subroutine source_virial_WL_MHD(sou_Tkin,sou_Pint,sou_Memf,sou_Wgra)
      real(8), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), &
      &                   sou_Memf(:,:,:), sou_Wgra(:,:,:)
    end subroutine source_virial_WL_MHD
  end interface
end module interface_source_virial_WL_MHD
