module interface_source_HaC_WL_EMF
  implicit none
  interface 
    subroutine source_HaC_WL_EMF(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_HaC_WL_EMF
  end interface
end module interface_source_HaC_WL_EMF
