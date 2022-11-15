subroutine interpo_gr2fl(grv,flv)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use interface_interpo_gr2fl_type0
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer  :: grv(:,:,:), flv(:,:,:)
  integer    :: irf, itf, ipf
  real(long) :: flvf
!
  flv(0:nrf,0:ntf,0:npf) = 0.0d0
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        call interpo_gr2fl_type0(flvf,grv,irf,itf,ipf)
        flv(irf,itf,ipf) = flvf
      end do
    end do
  end do
!
end subroutine interpo_gr2fl
