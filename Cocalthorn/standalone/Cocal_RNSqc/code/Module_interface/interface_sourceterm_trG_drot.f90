module interface_sourceterm_trG_drot
  implicit none
  interface 
    subroutine sourceterm_trG_drot(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_drot
  end interface
end module interface_sourceterm_trG_drot
