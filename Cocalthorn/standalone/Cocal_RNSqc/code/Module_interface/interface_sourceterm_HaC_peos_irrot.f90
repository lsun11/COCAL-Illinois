module interface_sourceterm_HaC_peos_irrot
  implicit none
  interface 
    subroutine sourceterm_HaC_peos_irrot(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_peos_irrot
  end interface
end module interface_sourceterm_HaC_peos_irrot
