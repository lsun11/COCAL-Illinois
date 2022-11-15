module interface_sourceterm_HaC_CF_peos_corot
  implicit none
  interface 
    subroutine sourceterm_HaC_CF_peos_corot(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_CF_peos_corot
  end interface
end module interface_sourceterm_HaC_CF_peos_corot
