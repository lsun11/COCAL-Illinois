module interface_hydrostatic_eq_peos_lecc
  implicit none
  interface 
    subroutine hydrostatic_eq_peos_lecc(emd)
      real(8), pointer :: emd(:,:,:)
    end subroutine hydrostatic_eq_peos_lecc
  end interface
end module interface_hydrostatic_eq_peos_lecc
