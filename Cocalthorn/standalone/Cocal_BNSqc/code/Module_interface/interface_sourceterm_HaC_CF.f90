module interface_sourceterm_HaC_CF
  implicit none
  interface 
    subroutine sourceterm_HaC_CF(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_CF
  end interface
end module interface_sourceterm_HaC_CF
