module interface_sourceterm_trG_WL
  implicit none
  interface 
    subroutine sourceterm_trG_WL(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_WL
  end interface
end module interface_sourceterm_trG_WL
