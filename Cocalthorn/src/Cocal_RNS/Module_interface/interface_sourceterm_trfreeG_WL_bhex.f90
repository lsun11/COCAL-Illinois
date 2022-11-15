module interface_sourceterm_trfreeG_WL_bhex
  implicit none
  interface 
    subroutine sourceterm_trfreeG_WL_bhex(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_WL_bhex
  end interface
end module interface_sourceterm_trfreeG_WL_bhex
