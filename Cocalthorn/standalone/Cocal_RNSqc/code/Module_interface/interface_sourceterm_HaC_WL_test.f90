module interface_sourceterm_HaC_WL_test
  implicit none
  interface 
    subroutine sourceterm_HaC_WL_test(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_WL_test
  end interface
end module interface_sourceterm_HaC_WL_test
