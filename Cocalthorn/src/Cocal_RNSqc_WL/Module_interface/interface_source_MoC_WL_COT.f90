module interface_source_MoC_WL_COT
  implicit none
  interface 
    subroutine source_MoC_WL_COT(souvec,cobj)
      real(8), pointer :: souvec(:,:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine source_MoC_WL_COT
  end interface
end module interface_source_MoC_WL_COT
