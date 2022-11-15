module interface_sourceterm_1bh_2pot_test
  implicit none
  interface 
    subroutine sourceterm_1bh_2pot_test(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_1bh_2pot_test
  end interface
end module interface_sourceterm_1bh_2pot_test
