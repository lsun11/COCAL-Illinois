module interface_source_trG_WL
  implicit none
  interface 
    subroutine source_trG_WL(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_trG_WL
  end interface
end module interface_source_trG_WL
