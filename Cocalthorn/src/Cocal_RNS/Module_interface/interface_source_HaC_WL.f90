module interface_source_HaC_WL
  implicit none
  interface 
    subroutine source_HaC_WL(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_HaC_WL
  end interface
end module interface_source_HaC_WL
