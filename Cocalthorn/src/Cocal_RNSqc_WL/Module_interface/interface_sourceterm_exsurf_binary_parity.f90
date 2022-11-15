module interface_sourceterm_exsurf_binary_parity
  implicit none
  interface 
    subroutine sourceterm_exsurf_binary_parity(fnc,sou_ex,dsou_ex,parchar)
      real(8), pointer :: fnc(:,:,:), sou_ex(:,:), dsou_ex(:,:)
      character(len=2) :: parchar
    end subroutine sourceterm_exsurf_binary_parity
  end interface
end module interface_sourceterm_exsurf_binary_parity
