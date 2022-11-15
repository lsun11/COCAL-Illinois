module interface_sourceterm_trfreeG_drot_SFC
  implicit none
  interface
    subroutine sourceterm_trfreeG_drot_SFC(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_drot_SFC
  end interface
end module interface_sourceterm_trfreeG_drot_SFC
