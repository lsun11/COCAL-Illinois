module interface_sourceterm_MoC_peos
  implicit none
  interface 
    subroutine sourceterm_MoC_peos(souvec,sou)
      real(8), pointer :: sou(:,:,:) 
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_peos
  end interface
end module interface_sourceterm_MoC_peos
