module interface_interpo_flsfc2flsph_midpoint
  implicit none
  interface 
    subroutine interpo_flsfc2flsph_midpoint(flv,flsphv)
      real(8), pointer :: flv(:,:,:), flsphv(:,:,:)
    end subroutine interpo_flsfc2flsph_midpoint
  end interface
end module interface_interpo_flsfc2flsph_midpoint
