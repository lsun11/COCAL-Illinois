subroutine vol_int_grav_bhex(sou,vol,irg_vol)
  use phys_constant, only  :   long
  use grid_parameter, only  :   nrg, ntg, npg
!  use weight_midpoint_grav, only : hwrtpg
  use grid_points_binary_excision,  only : ihpg_exin, ihpg_exout
  use weight_midpoint_binary_excision,  only : hwrtpg_ex

  implicit none
  integer     ::   ir,it,ip, ipin, ipout, irg_vol
  real(long)  ::   vol
  real(long), pointer  ::   sou(:,:,:)

  vol = 0.0d0
  do ir = 1, irg_vol
    do it = 1, ntg
      ipin = ihpg_exin(ir,it)
      ipout = ihpg_exout(ir,it)
      do ip = ipin, ipout
         vol = vol + sou(ir,it,ip)*hwrtpg_ex(ir,it,ip)
      end do
    end do
  end do
end subroutine vol_int_grav_bhex
