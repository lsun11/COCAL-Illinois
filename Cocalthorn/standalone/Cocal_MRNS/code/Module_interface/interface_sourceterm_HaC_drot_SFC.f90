module interface_sourceterm_HaC_drot_SFC
  implicit none
  interface 
    subroutine sourceterm_HaC_drot_SFC(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_drot_SFC
  end interface
end module interface_sourceterm_HaC_drot_SFC
