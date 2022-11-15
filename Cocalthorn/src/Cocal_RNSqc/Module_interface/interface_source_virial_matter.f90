module interface_source_virial_matter
  implicit none
  interface 
    subroutine source_virial_matter(sou_Tkin,sou_Pint)
      real(8), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:)
    end subroutine source_virial_matter
  end interface
end module interface_source_virial_matter
