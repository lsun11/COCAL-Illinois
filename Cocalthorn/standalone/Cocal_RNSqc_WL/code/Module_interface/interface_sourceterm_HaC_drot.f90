module interface_sourceterm_HaC_drot
  implicit none
  interface 
    subroutine sourceterm_HaC_drot(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC_drot
  end interface
end module interface_sourceterm_HaC_drot
