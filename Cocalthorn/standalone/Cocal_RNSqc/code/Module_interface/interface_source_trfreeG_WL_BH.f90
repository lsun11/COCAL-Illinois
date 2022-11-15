module interface_source_trfreeG_WL_BH
  implicit none
  interface
    subroutine source_trfreeG_WL_BH(souten)
      real(8), pointer :: souten(:,:,:,:)
    end subroutine source_trfreeG_WL_BH
  end interface
end module interface_source_trfreeG_WL_BH
