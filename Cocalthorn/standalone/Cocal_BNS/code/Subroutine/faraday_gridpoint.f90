subroutine faraday_gridpoint
  use phys_constant,  only : long
  use make_array_3d
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,  only : alph, bvxu, bvyu, bvzu
  use def_emfield, only : alva, vaxd, vayd, vazd
  use def_faraday_tensor, only : fxd, fyd, fzd, fxd_grid, fyd_grid, fzd_grid, &
  &                              fijd_grid, Lie_bFxd, Lie_bFyd, Lie_bFzd
  use def_emfield_derivatives, only : Lie_bAxd_grid, Lie_bAyd_grid, &
  &                                   Lie_bAzd_grid
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_midpoint_type0
  use interface_interpo_linear_type0
  implicit none
  integer :: irg, itg, ipg
  real(long) :: alpgp, ainvh, bvxugp, bvyugp, bvzugp
  real(long) :: vaxdgp, vaydgp, vazdgp
  real(long) :: fxdgc, fydgc, fzdgc, bvxugc, bvyugc, bvzugc
  real(long) :: dbvxdxgp, dbvxdygp, dbvxdzgp, &
  &             dbvydxgp, dbvydygp, dbvydzgp, &
  &             dbvzdxgp, dbvzdygp, dbvzdzgp
  real(long) :: dvaxdxgp, dvaxdygp, dvaxdzgp, &
  &             dvaydxgp, dvaydygp, dvaydzgp, &
  &             dvazdxgp, dvazdygp, dvazdzgp
  real(long) :: dbvxdxgc, dbvxdygc, dbvxdzgc, &
  &             dbvydxgc, dbvydygc, dbvydzgc, &
  &             dbvzdxgc, dbvzdygc, dbvzdzgc
  real(long) :: dfxdx, dfxdy, dfxdz, &
  &             dfydx, dfydy, dfydz, &
  &             dfzdx, dfzdy, dfzdz
  real(long) :: dalvadxgp, dalvadygp, dalvadzgp
  real(long) :: lie_bAx, lie_bAy, lie_bAz, lie_bFx, lie_bFy, lie_bFz
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_4th_gridpoint(bvxu,dbvxdxgp,dbvxdygp,dbvxdzgp,irg,itg,ipg)
        call grgrad_4th_gridpoint(bvyu,dbvydxgp,dbvydygp,dbvydzgp,irg,itg,ipg)
        call grgrad_4th_gridpoint(bvzu,dbvzdxgp,dbvzdygp,dbvzdzgp,irg,itg,ipg)
        call grgrad_4th_gridpoint(vaxd,dvaxdxgp,dvaxdygp,dvaxdzgp,irg,itg,ipg)
        call grgrad_4th_gridpoint(vayd,dvaydxgp,dvaydygp,dvaydzgp,irg,itg,ipg)
        call grgrad_4th_gridpoint(vazd,dvazdxgp,dvazdygp,dvazdzgp,irg,itg,ipg)
        call grgrad_4th_gridpoint(alva,dalvadxgp,dalvadygp,dalvadzgp &
        &                                                        ,irg,itg,ipg)
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
        Lie_bAxd_grid(irg,itg,ipg) = lie_bAx
        Lie_bAyd_grid(irg,itg,ipg) = lie_bAy
        Lie_bAzd_grid(irg,itg,ipg) = lie_bAz
        fxd_grid(irg,itg,ipg) = ainvh*(lie_bAx - dalvadxgp)
        fyd_grid(irg,itg,ipg) = ainvh*(lie_bAy - dalvadygp)
        fzd_grid(irg,itg,ipg) = ainvh*(lie_bAz - dalvadzgp)
        fijd_grid(irg,itg,ipg,1) = dvaydxgp - dvaxdygp
        fijd_grid(irg,itg,ipg,2) = dvazdxgp - dvaxdzgp
        fijd_grid(irg,itg,ipg,3) = dvazdygp - dvaydzgp
!
      end do
    end do
  end do
!
! compute Lie_bF at mid points
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call grgrad_midpoint_type0(bvxu,dbvxdxgc,dbvxdygc,dbvxdzgc,irg,itg,ipg)
        call grgrad_midpoint_type0(bvyu,dbvydxgc,dbvydygc,dbvydzgc,irg,itg,ipg)
        call grgrad_midpoint_type0(bvzu,dbvzdxgc,dbvzdygc,dbvzdzgc,irg,itg,ipg)
        call grgrad_midpoint_type0(fxd_grid,dfxdx,dfxdy,dfxdz,irg,itg,ipg)
        call grgrad_midpoint_type0(fyd_grid,dfydx,dfydy,dfydz,irg,itg,ipg)
        call grgrad_midpoint_type0(fzd_grid,dfzdx,dfzdy,dfzdz,irg,itg,ipg)
        call interpo_linear_type0(bvxugc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvyugc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzugc,bvzu,irg,itg,ipg)
        fxdgc = fxd(irg,itg,ipg)
        fydgc = fyd(irg,itg,ipg)
        fzdgc = fzd(irg,itg,ipg)
        lie_bFx = bvxugc*dfxdx   + bvyugc*dfxdy   + bvzugc*dfxdz &
        &       + fxdgc*dbvxdxgc + fydgc*dbvydxgc + fzdgc*dbvzdxgc
        lie_bFy = bvxugc*dfydx   + bvyugc*dfydy   + bvzugc*dfydz &
        &       + fxdgc*dbvxdygc + fydgc*dbvydygc + fzdgc*dbvzdygc
        lie_bFz = bvxugc*dfzdx   + bvyugc*dfzdy   + bvzugc*dfzdz &
        &       + fxdgc*dbvxdzgc + fydgc*dbvydzgc + fzdgc*dbvzdzgc
!
        Lie_bFxd(irg,itg,ipg) = lie_bFx
        Lie_bFyd(irg,itg,ipg) = lie_bFy
        Lie_bFzd(irg,itg,ipg) = lie_bFz
!
      end do
    end do
  end do
!
end subroutine faraday_gridpoint
