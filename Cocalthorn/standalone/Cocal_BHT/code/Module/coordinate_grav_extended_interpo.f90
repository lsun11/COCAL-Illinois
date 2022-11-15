! extended coordinate for the field
!______________________________________________
module coordinate_grav_extended_interpo
  use phys_constant, only : nnrg, nntg, nnpg, long
  implicit none
  real(long) ::   rgex_itp(-2:nnrg+2), &
  &              thgex_itp(-2:nntg+2), &
  &             phigex_itp(-2:nnpg+2)
  integer :: irgex_r_itp(-2:nnrg+2), &
  &          itgex_r_itp(0:nntg,-2:nnrg+2), &
  &          ipgex_r_itp(0:nnpg,-2:nnrg+2)
  integer :: itgex_th_itp(-2:nntg+2), &
  &          ipgex_th_itp(0:nnpg,-2:nntg+2)
  integer :: ipgex_phi_itp(-2:nnpg+2)
!
  real(long) :: hrgex_itp(-2:nnrg+2), &
  &            hthgex_itp(-2:nntg+2), &
  &           hphigex_itp(-2:nnpg+2)
  integer :: irgex_hr_itp(-2:nnrg+2), &
  &          itgex_hr_itp(1:nntg,-2:nnrg+2), &
  &          ipgex_hr_itp(1:nnpg,-2:nnrg+2)
  integer :: itgex_hth_itp(-2:nntg+2), &
  &          ipgex_hth_itp(1:nnpg,-2:nntg+2)
  integer :: ipgex_hphi_itp(-2:nnpg+2)
end module coordinate_grav_extended_interpo
