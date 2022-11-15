module interface_sourceterm_trG_WL_SEM
  implicit none
  interface 
    subroutine sourceterm_trG_WL_SEM(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG_WL_SEM
  end interface
end module interface_sourceterm_trG_WL_SEM
