module interface_sourceterm_outsurf_interpo_from_asympto_parity_mpt
  implicit none
  interface 
    subroutine sourceterm_outsurf_interpo_from_asympto_parity_mpt &
    &              (impt_bin,impt_apt,fnc,sou_out,dsou_out,parchar)
      real(8), pointer :: fnc(:,:,:), sou_out(:,:), dsou_out(:,:)
      integer :: impt_bin, impt_apt
      character(len=2) :: parchar
    end subroutine sourceterm_outsurf_interpo_from_asympto_parity_mpt
  end interface
end module interface_sourceterm_outsurf_interpo_from_asympto_parity_mpt
