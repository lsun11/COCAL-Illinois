module interface_sourceterm_insurf_asymptotic_patch
  implicit none
  interface 
    subroutine sourceterm_insurf_asymptotic_patch(impt_bin,impt_apt,fnc, &
    &                                             sou_in,dsou_in)
      real(8), pointer :: fnc(:,:,:), sou_in(:,:), dsou_in(:,:)
      integer, intent(IN) :: impt_bin, impt_apt
    end subroutine sourceterm_insurf_asymptotic_patch
  end interface
end module interface_sourceterm_insurf_asymptotic_patch
