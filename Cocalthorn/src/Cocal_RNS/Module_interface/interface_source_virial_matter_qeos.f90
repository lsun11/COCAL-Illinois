module interface_source_virial_matter_qeos
  implicit none
  interface 
    subroutine source_virial_matter_qeos(sou_Tkin,sou_Pint)
      real(8), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:)
    end subroutine source_virial_matter_qeos
  end interface
end module interface_source_virial_matter_qeos
