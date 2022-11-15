subroutine interpo_gr2fl_surface(grv,flv)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use interface_interpo_gr2fl_type0
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer  :: grv(:,:,:), flv(:,:)
  integer    :: irf, itf, ipf
  real(long) :: flvf
!
  flv(0:ntf,0:npf) = 0.0d0
!
  irf = nrf
  do ipf = 0, npf
    do itf = 0, ntf
      call interpo_gr2fl_type0(flvf,grv,irf,itf,ipf)
      flv(itf,ipf) = flvf
    end do
  end do
!
end subroutine interpo_gr2fl_surface
