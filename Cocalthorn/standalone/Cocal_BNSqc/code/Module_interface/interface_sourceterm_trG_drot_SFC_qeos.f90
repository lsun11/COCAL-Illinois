module interface_sourceterm_trG_drot_SFC_qeos
  implicit none
  interface 
    subroutine sourceterm_trG_drot_SFC_qeos(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_drot_SFC_qeos
  end interface
end module interface_sourceterm_trG_drot_SFC_qeos
