module interface_sourceterm_trfreeG_WL_SEM
  implicit none
  interface 
    subroutine sourceterm_trfreeG_WL_SEM(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_WL_SEM
  end interface
end module interface_sourceterm_trfreeG_WL_SEM
