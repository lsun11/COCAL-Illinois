subroutine calc_weight_midpoint_fluid_sphcoord
  use grid_parameter,       only : nrf, ntg, npg
  use coordinate_grav_r,    only : rg, hrg
  use weight_midpoint_grav, only : hwdrg, hwdtg, hwdpg
  use def_matter,           only : rs
  use weight_midpoint_fluid_sphcoord
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long) :: wgr, wgt, wgp, rrff
  integer :: irg, itg, ipg
!
! --  weight of integlation for volume integral in GR coordinate,
! --  with modified weight for the stellar surface
!
!  hwrtpg_fc(:,:,:) = 0.0d0
  hwrtpg_fc = 0.0d0

  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(rrff,rs,itg,ipg)
      do irg = 1, nrf
        wgr = hwdrg(irg)
        wgt = hwdtg(itg)
        wgp = hwdpg(ipg)
        hwrtpg_fc(irg,itg,ipg) = wgr*wgt*wgp
        if (hrg(irg+1).gt.rrff*rg(nrf)) then
          wgr = hrg(irg)**2*(rrff*rg(nrf) - rg(irg-1))
          hwrtpg_fc(irg,itg,ipg) = wgr*wgt*wgp
          exit
        end if
      end do
    end do
  end do
! 
!  do ipg = 1, npg
!    do itg = 1, ntg
!      write(6,*)  ipg,itg, hwrtpg_fc(nrf-1,itg,ipg), hwrtpg_fc(nrf,itg,ipg)
!    end do
!  end do

end subroutine calc_weight_midpoint_fluid_sphcoord

