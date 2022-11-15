module interface_sourceterm_AHfinder_WL
  implicit none
  interface 
    subroutine sourceterm_AHfinder_WL(sou)
      real(8), pointer :: sou(:,:)
    end subroutine sourceterm_AHfinder_WL
  end interface
end module interface_sourceterm_AHfinder_WL
