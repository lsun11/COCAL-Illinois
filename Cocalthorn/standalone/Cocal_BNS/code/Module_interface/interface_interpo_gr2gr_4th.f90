module interface_interpo_gr2gr_4th
  implicit none
  interface 
    subroutine interpo_gr2gr_4th(fnc,cfn,rc,thc,phic)
      real(8), pointer     :: fnc(:,:,:)
      real(8), intent(out) :: cfn
      real(8) ::  rc, thc, phic
    end subroutine interpo_gr2gr_4th
  end interface
end module interface_interpo_gr2gr_4th
