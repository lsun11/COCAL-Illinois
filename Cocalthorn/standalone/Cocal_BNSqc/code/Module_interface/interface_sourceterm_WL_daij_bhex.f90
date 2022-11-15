module interface_sourceterm_WL_daij_bhex
  implicit none
  interface 
    subroutine sourceterm_WL_daij_bhex(sou_baij)
      real(8), pointer :: sou_baij(:,:,:)
    end subroutine sourceterm_WL_daij_bhex
  end interface
end module interface_sourceterm_WL_daij_bhex
