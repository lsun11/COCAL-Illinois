subroutine update_matter(potf,mtfield,convf)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  implicit none
  real(long), pointer    :: potf(:,:,:)
  real(long), pointer    :: mtfield(:,:,:)
  real(long), intent(in) :: convf
  integer    :: irf, itf, ipf
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        mtfield(irf,itf,ipf) = convf *   potf(irf,itf,ipf) &
      &                + (1.d0-convf)*mtfield(irf,itf,ipf)
      end do
    end do
  end do
end subroutine update_matter
