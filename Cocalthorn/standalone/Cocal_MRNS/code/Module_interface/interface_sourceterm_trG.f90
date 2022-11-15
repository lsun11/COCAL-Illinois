module interface_sourceterm_trG
  implicit none
  interface 
    subroutine sourceterm_trG(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG
  end interface
end module interface_sourceterm_trG
