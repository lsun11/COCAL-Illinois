module interface_source_admmom_asympto
  implicit none
  interface 
    subroutine source_admmom_asympto(sousfv,irg)
      real(8), pointer :: sousfv(:,:,:)
      integer          :: irg
    end subroutine source_admmom_asympto
  end interface
end module interface_source_admmom_asympto
