module interface_sourceterm_HaC
  implicit none
  interface 
    subroutine sourceterm_HaC(sou)    
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC
  end interface
end module interface_sourceterm_HaC
