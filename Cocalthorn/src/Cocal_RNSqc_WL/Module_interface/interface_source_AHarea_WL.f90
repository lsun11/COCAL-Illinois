module interface_source_AHarea_WL
  implicit none
  interface 
    subroutine source_AHarea_WL(sou)
      real(8), pointer :: sou(:,:)
    end subroutine source_AHarea_WL
  end interface
end module interface_source_AHarea_WL
