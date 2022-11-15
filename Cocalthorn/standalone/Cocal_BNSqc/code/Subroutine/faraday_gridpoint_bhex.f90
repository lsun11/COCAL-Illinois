subroutine faraday_gridpoint
  use phys_constant,  only : long
  use make_array_3d
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,  only : alph, bvxu, bvyu, bvzu
  use def_emfield,  only : alva, vaxd, vayd, vazd
  use def_faraday_tensor,  only : fxd, fyd, fzd
  use interface_grgrad_4th_gridpoint_bhex
  implicit none
  integer :: irg, itg, ipg
  real(long) :: alpgp, ainvh, bvxugp, bvyugp, bvzugp
  real(long) :: vaxdgp, vaydgp, vazdgp
  real(long) :: dbvxdxgp, dbvxdygp, dbvxdzgp &
  &             dbvydxgp, dbvydygp, dbvydzgp &
  &             dbvzdxgp, dbvzdygp, dbvzdzgp
  real(long) :: dvaxdxgp, dvaxdygp, dvaxdzgp &
  &             dvaydxgp, dvaydygp, dvaydzgp &
  &             dvazdxgp, dvazdygp, dvazdzgp
  real(long) :: dalvadxgp, dalvadygp, dalvadzgp
  real(long) :: lie_bAx, lie_bAy, lie_bAz
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_4th_gridpoint_bhex(bvxu,dbvxdxgp,dbvxdygp,dbvxdzgp, &
        &                                                 ,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(bvyu,dbvydxgp,dbvydygp,dbvydzgp, &
        &                                                 ,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(bvzu,dbvzdxgp,dbvzdygp,dbvzdzgp, &
        &                                                 ,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(vaxd,dvaxdxgp,dvaxdygp,dvaxdzgp, &
        &                                                 ,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(vayd,dvaydxgp,dvaydygp,dvaydzgp, &
        &                                                 ,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(vazd,dvazdxgp,dvazdygp,dvazdzgp, &
        &                                                 ,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(alva,dalvadxgp,dalvadygp,dalvadzgp, &
        &                                                 ,irg,itg,ipg)
        ainvh = 1.0d0/alph(irg,itg,ipg)
        bvxugp = bvxu(irg,itg,ipg)
        bvyugp = bvyu(irg,itg,ipg)
        bvzugp = bvzu(irg,itg,ipg)
        vaxdgp = vaxd(irg,itg,ipg)
        vaydgp = vayd(irg,itg,ipg)
        vazdgp = vazd(irg,itg,ipg)
!
        lie_bAx = bvxugp*dvaxdxgp + bvyugp*dvaxdygp + bvzugp*dvaxdzgp &
        &       + vaxdgp*dbvxdxgp + vaydgp*dbvydxgp + vazdgp*dbvzdxgp
        lie_bAy = bvxugp*dvaydxgp + bvyugp*dvaydygp + bvzugp*dvaydzgp &
        &       + vaxdgp*dbvxdygp + vaydgp*dbvydygp + vazdgp*dbvzdygp
        lie_bAz = bvxugp*dvazdxgp + bvyugp*dvazdygp + bvzugp*dvazdzgp &
        &       + vaxdgp*dbvxdzgp + vaydgp*dbvydzgp + vazdgp*dbvzdzgp
!
        fxd_grid(irg,itg,ipg) = ainvh*(lie_bAx - dalvadxgp)
        fyd_grid(irg,itg,ipg) = ainvh*(lie_bAy - dalvadygp)
        fzd_grid(irg,itg,ipg) = ainvh*(lie_bAz - dalvadzgp)
!
      end do
    end do
  end do
!
end subroutine faraday_gridpoint
