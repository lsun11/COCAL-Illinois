subroutine update_matter(potf,mtfield,convf)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  implicit none
  real(long), pointer    :: potf(:,:,:)
  real(long), pointer    :: mtfield(:,:,:)
  real(long), intent(in) :: convf
  real(long)             :: rand, frac, convrand
  integer                :: irf,itf,ipf
  mtfield(0:nrf,0:ntf,0:npf) =       convf *   potf(0:nrf,0:ntf,0:npf) &
  &                          + (1.d0-convf)*mtfield(0:nrf,0:ntf,0:npf)
!  mtfield(0:nrf-2,0:ntf,0:npf) =       convf *   potf(0:nrf-2,0:ntf,0:npf) &
!  &                            + (1.d0-convf)*mtfield(0:nrf-2,0:ntf,0:npf)
!no  do ipf = 0, npf
!no    do itf = 0, ntf
!no!      do irf = nrf-1, nrf
!no      do irf = 0, nrf
!no        call random_number(rand)
!no!no        rand  = rand - 0.5d0
!no!no        frac = 0.25d0*convf
!no!no!        convrand = convf - 0.5*frac + frac*rand
!no!no        convrand = convf + frac*rand
!no!no        mtfield(irf,itf,ipf) = (     convrand)*   potf(irf,itf,ipf) &
!no!no        &                    + (1.d0-convrand)*mtfield(irf,itf,ipf)
!no        if (rand.le.0.5) then
!no          mtfield(irf,itf,ipf) =       convf *   potf(irf,itf,ipf) &
!no          &                    + (1.d0-convf)*mtfield(irf,itf,ipf)
!no        end if
!no      end do
!no    end do
!no  end do
end subroutine update_matter
