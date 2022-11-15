module interface_source_HaC_CF
  implicit none
  interface 
    subroutine source_HaC_CF(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_HaC_CF
  end interface
end module interface_source_HaC_CF
