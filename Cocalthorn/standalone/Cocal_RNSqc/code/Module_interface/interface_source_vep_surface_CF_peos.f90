module interface_source_vep_surface_CF_peos
  implicit none
  interface 
    subroutine source_vep_surface_CF_peos(vpot_v,surp)
      real(8), pointer  ::   surp(:,:), vpot_v(:,:,:)
    end subroutine source_vep_surface_CF_peos
  end interface
end module interface_source_vep_surface_CF_peos
