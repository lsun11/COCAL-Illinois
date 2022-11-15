module interface_source_trfreeG_WL_EMF
  implicit none
  interface
    subroutine source_trfreeG_WL_EMF(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine source_trfreeG_WL_EMF
  end interface
end module interface_source_trfreeG_WL_EMF
