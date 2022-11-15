module interface_hydro_irbns_vep_CF_peos_lecc
  implicit none
  interface 
    subroutine hydro_irbns_vep_CF_peos_lecc(vpot,iter,impt,ihy)
      real(8), pointer :: vpot(:,:,:)
      integer :: iter,impt,ihy
    end subroutine hydro_irbns_vep_CF_peos_lecc
  end interface
end module interface_hydro_irbns_vep_CF_peos_lecc
