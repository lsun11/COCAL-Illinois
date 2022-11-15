module interface_sourceterm_trfreeG_WL
  implicit none
  interface 
    subroutine sourceterm_trfreeG_WL(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_WL
  end interface
end module interface_sourceterm_trfreeG_WL
