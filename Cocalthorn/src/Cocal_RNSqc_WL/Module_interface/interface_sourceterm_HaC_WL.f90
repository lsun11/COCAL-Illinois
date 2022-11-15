module interface_sourceterm_HaC_WL
  implicit none
  interface 
    subroutine sourceterm_HaC_WL(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_WL
  end interface
end module interface_sourceterm_HaC_WL
