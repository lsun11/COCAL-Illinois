module interface_sourceterm_insurf_asympto_interpo_from_parity_mpt
  implicit none
  interface 
    subroutine sourceterm_insurf_asympto_interpo_from_parity_mpt &
    &                        (impt_bin,impt_apt,fnc,sou_in,dsou_in,parchar)
      real(8), pointer :: fnc(:,:,:), sou_in(:,:), dsou_in(:,:)
      integer, intent(in) :: impt_bin, impt_apt
      character(len=2) :: parchar
    end subroutine sourceterm_insurf_asympto_interpo_from_parity_mpt
  end interface
end module interface_sourceterm_insurf_asympto_interpo_from_parity_mpt
