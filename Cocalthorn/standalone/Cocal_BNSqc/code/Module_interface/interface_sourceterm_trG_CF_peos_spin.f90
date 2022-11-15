module interface_sourceterm_trG_CF_peos_spin
  implicit none
  interface 
    subroutine sourceterm_trG_CF_peos_spin(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_CF_peos_spin
  end interface
end module interface_sourceterm_trG_CF_peos_spin
