module interface_source_HaC_CF_qeos
  implicit none
  interface 
    subroutine source_HaC_CF_qeos(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_HaC_CF_qeos
  end interface
end module interface_source_HaC_CF_qeos
