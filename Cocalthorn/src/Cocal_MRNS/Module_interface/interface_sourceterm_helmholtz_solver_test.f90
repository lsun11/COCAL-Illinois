module interface_sourceterm_helmholtz_solver_test
  implicit none
  interface 
    subroutine sourceterm_helmholtz_solver_test(sou,char_sou)
      real(8), pointer :: sou(:,:,:)
      character(2)     :: char_sou
    end subroutine sourceterm_helmholtz_solver_test
  end interface
end module interface_sourceterm_helmholtz_solver_test
