module interface_source_vep_CF_peos_lecc
  implicit none
  interface 
    subroutine source_vep_CF_peos_lecc(souv)
      real(8), pointer :: souv(:,:,:)
    end subroutine source_vep_CF_peos_lecc
  end interface
end module interface_source_vep_CF_peos_lecc
