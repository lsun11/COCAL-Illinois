module interface_source_trfreeG_WL
  implicit none
  interface
    subroutine source_trfreeG_WL(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine source_trfreeG_WL
  end interface
end module interface_source_trfreeG_WL
