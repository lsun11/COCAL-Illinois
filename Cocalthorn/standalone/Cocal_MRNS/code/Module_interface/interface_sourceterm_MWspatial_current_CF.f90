module interface_sourceterm_MWspatial_current_CF
  implicit none
  interface 
    subroutine sourceterm_MWspatial_current_CF(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MWspatial_current_CF
  end interface
end module interface_sourceterm_MWspatial_current_CF
