module interface_sourceterm_trG_peos
  implicit none
  interface 
    subroutine sourceterm_trG_peos(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_peos
  end interface
end module interface_sourceterm_trG_peos
