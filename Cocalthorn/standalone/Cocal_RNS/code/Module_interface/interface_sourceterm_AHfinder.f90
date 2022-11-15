module interface_sourceterm_AHfinder
  implicit none
  interface 
    subroutine sourceterm_AHfinder(sou)
      real(8), pointer :: sou(:,:)
    end subroutine sourceterm_AHfinder
  end interface
end module interface_sourceterm_AHfinder
