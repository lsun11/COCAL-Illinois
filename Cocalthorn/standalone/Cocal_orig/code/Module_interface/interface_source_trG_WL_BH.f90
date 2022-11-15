module interface_source_trG_WL_BH
  implicit none
  interface 
    subroutine source_trG_WL_BH(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_trG_WL_BH
  end interface
end module interface_source_trG_WL_BH
