module interface_interpo_gr2fl_type0
  implicit none
  interface 
    subroutine interpo_gr2fl_type0(flv,grv,irf,itf,ipf)
      real(8) :: flv
      real(8), pointer :: grv(:,:,:)
      integer :: irf, itf, ipf
    end subroutine interpo_gr2fl_type0
  end interface
end module interface_interpo_gr2fl_type0
