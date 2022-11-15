subroutine vol_int_fluid(souf,vol)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrf, ntf, npf
  use coordinate_grav_r, only  :   drg, rg
  use coordinate_grav_theta, only  : dthg
  use coordinate_grav_phi, only  :   dphig
  use trigonometry_grav_theta, only  : sinthg, hsinthg
  use weight_midpoint_fluid, only : hwrtpf, siwrtpf
  use def_matter, only : rs
!  use interface_interpo_linear_type0
  use interface_interpo_linear_type0_2Dsurf
  implicit none
!
  integer     ::   ir,it,ip
  real(long)  ::   vol, hsouf, hsurf
  real(long),pointer  ::   souf(:,:,:)
!
! -- Mid-point integration overestimates the mass by 1% for nr = 16.
!  vol = 0.0d0
!  do ip = 1, npf
!    do it = 1, ntf
!      do ir = 1, nrf
!!!      hsouf = 0.125*(souf(ir,it,ip) + souf(ir-1,it,ip) + souf(ir,it-1,ip) &
!!!            + souf(ir-1,it-1,ip) + souf(ir,it,ip-1) + souf(ir-1,it,ip-1) & 
!!!            + souf(ir,it-1,ip-1) + souf(ir-1,it-1,ip-1))
!        call interpo_linear_type0(hsouf,souf,ir,it,ip)
!        call interpo_linear_type0_2Dsurf(hsurf,rs,it,ip)
!        vol = vol + hsouf * hsurf**3*hwrtpg(ir,it,ip)
!      end do
!    end do
!  end do
!
  vol = 0.0d0
  do ip = 1, npf
    do it = 1, ntf
      do ir = 0, nrf
        hsouf = 0.25*(souf(ir,it,ip)   + souf(ir,it-1,ip) &
                    + souf(ir,it,ip-1) + souf(ir,it-1,ip-1))
        call interpo_linear_type0_2Dsurf(hsurf,rs,it,ip)
        vol = vol + hsouf * hsurf**3 *siwrtpf(ir,it,ip)
      end do
    end do
  end do
!
end subroutine vol_int_fluid
