module interface_source_AHarea_CF
  implicit none
  interface 
    subroutine source_AHarea_CF(sou)
      real(8), pointer :: sou(:,:)
    end subroutine source_AHarea_CF
  end interface
end module interface_source_AHarea_CF
