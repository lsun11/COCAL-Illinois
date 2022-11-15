module interface_sourceterm_outsurf_interpo_from_asympto_mpt
  implicit none
  interface 
    subroutine sourceterm_outsurf_interpo_from_asympto_mpt &
    &              (impt_bin,impt_apt,fnc,sou_out,dsou_out)
      real(8), pointer :: fnc(:,:,:), sou_out(:,:), dsou_out(:,:)
      integer :: impt_bin, impt_apt
    end subroutine sourceterm_outsurf_interpo_from_asympto_mpt
  end interface
end module interface_sourceterm_outsurf_interpo_from_asympto_mpt
