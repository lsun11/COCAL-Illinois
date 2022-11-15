subroutine vol_int_fluid(souf,vol)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrf, ntf, npf
  use coordinate_grav_r, only  :   drg, rg
  use coordinate_grav_theta, only  : dthg
  use coordinate_grav_phi, only  :   dphig
  use trigonometry_grav_theta, only  : sinthg, hsinthg
  use weight_midpoint_fluid, only : rtsiwrtpf
  use def_matter, only : rs
  implicit none
!
  integer     ::   ir,it,ip
  real(long)  ::   vol, hsouf, hsurf
  real(long),pointer  ::   souf(:,:,:)
!
  vol = 0.0d0
  do ip = 0, npf
    do it = 0, ntf
      hsurf = rs(it,ip)
      do ir = 0, nrf
        hsouf = souf(ir,it,ip)
        vol = vol + hsouf * hsurf**3 *rtsiwrtpf(ir,it,ip)
      end do
    end do
  end do
!
end subroutine vol_int_fluid
