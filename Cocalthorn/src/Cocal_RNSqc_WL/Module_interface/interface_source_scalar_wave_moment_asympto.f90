module interface_source_scalar_wave_moment_asympto
  implicit none
  interface 
    subroutine source_scalar_wave_moment_asympto(sousfv,irg)
      real(8), pointer :: sousfv(:,:,:)
      integer          :: irg
    end subroutine source_scalar_wave_moment_asympto
  end interface
end module interface_source_scalar_wave_moment_asympto
