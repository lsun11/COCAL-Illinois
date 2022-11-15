module interface_source_MWtemp_WL
  implicit none
  interface 
    subroutine source_MWtemp_WL(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_MWtemp_WL
  end interface
end module interface_source_MWtemp_WL
