module interface_interpo_flsph2flsfc
  implicit none
  interface 
    subroutine interpo_flsph2flsfc(flsphv,flv)
      real(8), pointer :: flv(:,:,:), flsphv(:,:,:)
    end subroutine interpo_flsph2flsfc
  end interface
end module interface_interpo_flsph2flsfc
