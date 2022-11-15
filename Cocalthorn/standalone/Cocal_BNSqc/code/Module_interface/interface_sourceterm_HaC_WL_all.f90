module interface_sourceterm_HaC_WL_all
  implicit none
  interface 
    subroutine sourceterm_HaC_WL_all(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_WL_all
  end interface
end module interface_sourceterm_HaC_WL_all
