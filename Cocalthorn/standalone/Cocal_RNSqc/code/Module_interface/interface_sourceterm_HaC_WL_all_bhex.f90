module interface_sourceterm_HaC_WL_all_bhex
  implicit none
  interface 
    subroutine sourceterm_HaC_WL_all_bhex(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_WL_all_bhex
  end interface
end module interface_sourceterm_HaC_WL_all_bhex
