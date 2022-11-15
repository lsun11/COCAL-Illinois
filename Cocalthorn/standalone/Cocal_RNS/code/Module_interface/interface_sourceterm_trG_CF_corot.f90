module interface_sourceterm_trG_CF_corot
  implicit none
  interface 
    subroutine sourceterm_trG_CF_corot(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_CF_corot
  end interface
end module interface_sourceterm_trG_CF_corot
