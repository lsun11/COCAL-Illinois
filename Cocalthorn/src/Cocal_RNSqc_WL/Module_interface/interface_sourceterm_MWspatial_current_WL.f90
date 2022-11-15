module interface_sourceterm_MWspatial_current_WL
  implicit none
  interface 
    subroutine sourceterm_MWspatial_current_WL(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MWspatial_current_WL
  end interface
end module interface_sourceterm_MWspatial_current_WL
