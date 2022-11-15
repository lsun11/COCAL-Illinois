module interface_sourceterm_qij_WL_bhex
  implicit none
  interface 
    subroutine sourceterm_qij_WL_bhex(iter_count,sou_qij)
      real(8), pointer :: sou_qij(:,:,:,:)
      integer :: iter_count
    end subroutine sourceterm_qij_WL_bhex
  end interface
end module interface_sourceterm_qij_WL_bhex
