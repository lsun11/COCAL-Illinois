module interface_source_fnc_moment_asympto
  implicit none
  interface 
    subroutine source_fnc_moment_asympto(fnc,sousfv,irg)
      real(8), pointer :: fnc(:,:,:), sousfv(:,:,:)
      integer          :: irg
    end subroutine source_fnc_moment_asympto
  end interface
end module interface_source_fnc_moment_asympto
