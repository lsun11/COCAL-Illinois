module interface_source_trG_CF
  implicit none
  interface 
    subroutine source_trG_CF(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_trG_CF
  end interface
end module interface_source_trG_CF
