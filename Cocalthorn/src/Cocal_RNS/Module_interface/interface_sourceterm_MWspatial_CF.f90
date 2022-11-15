module interface_sourceterm_MWspatial_CF
  implicit none
  interface 
    subroutine sourceterm_MWspatial_CF(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MWspatial_CF
  end interface
end module interface_sourceterm_MWspatial_CF
