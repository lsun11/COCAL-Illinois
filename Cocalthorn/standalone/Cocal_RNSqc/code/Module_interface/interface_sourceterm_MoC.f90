module interface_sourceterm_MoC
  implicit none
  interface 
    subroutine sourceterm_MoC(souvec,sou)
      real(8), pointer :: sou(:,:,:) 
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC
  end interface
end module interface_sourceterm_MoC
