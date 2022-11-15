module interface_sourceterm_insurf_asympto_interpo_from_mpt
  implicit none
  interface 
    subroutine sourceterm_insurf_asympto_interpo_from_mpt &
    &                        (impt_bin,impt_apt,fnc,sou_in,dsou_in)
      real(8), pointer :: fnc(:,:,:), sou_in(:,:), dsou_in(:,:)
      integer, intent(in) :: impt_bin, impt_apt
    end subroutine sourceterm_insurf_asympto_interpo_from_mpt
  end interface
end module interface_sourceterm_insurf_asympto_interpo_from_mpt
