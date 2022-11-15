module interface_sourceterm_MWtemp_CF
  implicit none
  interface 
    subroutine sourceterm_MWtemp_CF(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_MWtemp_CF
  end interface
end module interface_sourceterm_MWtemp_CF
