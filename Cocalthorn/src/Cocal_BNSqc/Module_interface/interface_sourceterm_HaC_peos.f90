module interface_sourceterm_HaC_peos
  implicit none
  interface 
    subroutine sourceterm_HaC_peos(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_peos
  end interface
end module interface_sourceterm_HaC_peos
