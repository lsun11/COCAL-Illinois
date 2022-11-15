module interface_source_trfreeG_WL_COT
  implicit none
  interface
    subroutine source_trfreeG_WL_COT(souten,cobj)
      real(8), pointer :: souten(:,:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine source_trfreeG_WL_COT
  end interface
end module interface_source_trfreeG_WL_COT
