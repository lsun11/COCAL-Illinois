module interface_sourceterm_trfreeG_WL_EMF
  implicit none
  interface 
    subroutine sourceterm_trfreeG_WL_EMF(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_WL_EMF
  end interface
end module interface_sourceterm_trfreeG_WL_EMF
