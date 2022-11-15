module interface_source_trG_WL_COT
  implicit none
  interface 
    subroutine source_trG_WL_COT(sou,cobj)
      real(8), pointer :: sou(:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine source_trG_WL_COT
  end interface
end module interface_source_trG_WL_COT
