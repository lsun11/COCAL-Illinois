module interface_sourceterm_trfreeG_corot
  implicit none
  interface
    subroutine sourceterm_trfreeG_corot(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine sourceterm_trfreeG_corot
  end interface
end module interface_sourceterm_trfreeG_corot
