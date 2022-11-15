module interface_sourceterm_trG_CF_pBH
  implicit none
  interface 
    subroutine sourceterm_trG_CF_pBH(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_CF_pBH
  end interface
end module interface_sourceterm_trG_CF_pBH
