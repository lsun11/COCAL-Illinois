module interface_sourceterm_MWtemp_current
  implicit none
  interface 
    subroutine sourceterm_MWtemp_current(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_MWtemp_current
  end interface
end module interface_sourceterm_MWtemp_current
