module interface_sourceterm_trG_WL_EMF
  implicit none
  interface 
    subroutine sourceterm_trG_WL_EMF(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_WL_EMF
  end interface
end module interface_sourceterm_trG_WL_EMF
