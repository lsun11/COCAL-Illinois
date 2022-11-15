module interface_sourceterm_trG_CF_helical
  implicit none
  interface 
    subroutine sourceterm_trG_CF_helical(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_CF_helical
  end interface
end module interface_sourceterm_trG_CF_helical
