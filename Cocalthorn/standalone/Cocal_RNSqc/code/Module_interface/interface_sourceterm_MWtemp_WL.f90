module interface_sourceterm_MWtemp_WL
  implicit none
  interface 
    subroutine sourceterm_MWtemp_WL(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_MWtemp_WL
  end interface
end module interface_sourceterm_MWtemp_WL
