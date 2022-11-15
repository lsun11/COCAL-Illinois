module interface_sourceterm_trfreeG_WL_BHT
  implicit none
  interface 
    subroutine sourceterm_trfreeG_WL_BHT(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_WL_BHT
  end interface
end module interface_sourceterm_trfreeG_WL_BHT
