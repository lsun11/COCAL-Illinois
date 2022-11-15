module interface_sourceterm_trG_CF_peos_irrot
  implicit none
  interface 
    subroutine sourceterm_trG_CF_peos_irrot(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_CF_peos_irrot
  end interface
end module interface_sourceterm_trG_CF_peos_irrot
