subroutine liegmab_gridpoint
  use grid_parameter, only : nrg, ntg, npg
  use def_Lie_derivatives_grid, only : elpxx_grid, elpxy_grid, elpxz_grid, &
  &                               elpyy_grid, elpyz_grid, elpzz_grid, &
  &                               rlpxx_grid, rlpxy_grid, rlpxz_grid, &
  &                               rlpyy_grid, rlpyz_grid, rlpzz_grid
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_dvphi, only : dphiu
  use interface_grdphi_gridpoint_type0
  implicit none
  real(8) :: dhxxdp, dhxxdx, dhxxdy, dhxxdz, &
  &          dhxydp, dhxydx, dhxydy, dhxydz, &
  &          dhxzdp, dhxzdx, dhxzdy, dhxzdz, &
  &          dhyxdp, dhyxdx, dhyxdy, dhyxdz, &
  &          dhyydp, dhyydx, dhyydy, dhyydz, &
  &          dhyzdp, dhyzdx, dhyzdy, dhyzdz, &
  &          dhzxdp, dhzxdx, dhzxdy, dhzxdz, &
  &          dhzydp, dhzydx, dhzydy, dhzydz, &
  &          dhzzdp, dhzzdx, dhzzdy, dhzzdz, &
  &          fa13, trlie, &
  &          gmxxd, gmxyd, gmxzd, gmxxu, gmxyu, gmxzu, &
  &          gmyxd, gmyyd, gmyzd, gmyxu, gmyyu, gmyzu, &
  &          gmzxd, gmzyd, gmzzd, gmzxu, gmzyu, gmzzu
  integer :: ipg, irg, itg
!
! --- Compute Lie derivative of tgamma w.r.t. phi, and 
! --- and conformally invariant derivative of tgamma w.r.t. phi. 
! NOTE: dphiu(a,b) = D_b phi^a
!
  fa13 = 1.0d0/3.0d0
  dphiu(1:3,1:3) = 0.0d0
  dphiu(1,2) =-1.0d0
  dphiu(2,1) = 1.0d0
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        gmxxd=hxxd(irg,itg,ipg)
        gmxyd=hxyd(irg,itg,ipg)
        gmxzd=hxzd(irg,itg,ipg)
        gmyyd=hyyd(irg,itg,ipg)
        gmyzd=hyzd(irg,itg,ipg)
        gmzzd=hzzd(irg,itg,ipg)
        gmxxd = gmxxd + 1.0d0
        gmyyd = gmyyd + 1.0d0
        gmzzd = gmzzd + 1.0d0
        gmyxd = gmxyd
        gmzxd = gmxzd
        gmzyd = gmyzd
        gmxxu=hxxu(irg,itg,ipg)
        gmxyu=hxyu(irg,itg,ipg)
        gmxzu=hxzu(irg,itg,ipg)
        gmyyu=hyyu(irg,itg,ipg)
        gmyzu=hyzu(irg,itg,ipg)
        gmzzu=hzzu(irg,itg,ipg)
        gmxxu = gmxxu + 1.0d0
        gmyyu = gmyyu + 1.0d0
        gmzzu = gmzzu + 1.0d0
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu
!
        call grdphi_gridpoint_type0(hxxd,dhxxdp,irg,itg,ipg)
        call grdphi_gridpoint_type0(hxyd,dhxydp,irg,itg,ipg)
        call grdphi_gridpoint_type0(hxzd,dhxzdp,irg,itg,ipg)
        call grdphi_gridpoint_type0(hyyd,dhyydp,irg,itg,ipg)
        call grdphi_gridpoint_type0(hyzd,dhyzdp,irg,itg,ipg)
        call grdphi_gridpoint_type0(hzzd,dhzzdp,irg,itg,ipg)
!
        rlpxx_grid(irg,itg,ipg) = dhxxdp + gmxyd*dphiu(2,1) + gmyxd*dphiu(2,1)
        rlpxy_grid(irg,itg,ipg) = dhxydp + gmxxd*dphiu(1,2) + gmyyd*dphiu(2,1)
        rlpxz_grid(irg,itg,ipg) = dhxzdp + gmyzd*dphiu(2,1)
        rlpyy_grid(irg,itg,ipg) = dhyydp + gmyxd*dphiu(1,2) + gmxyd*dphiu(1,2)
        rlpyz_grid(irg,itg,ipg) = dhyzdp + gmxzd*dphiu(1,2)
        rlpzz_grid(irg,itg,ipg) = dhzzdp
!
        dhxxdp = rlpxx_grid(irg,itg,ipg)
        dhxydp = rlpxy_grid(irg,itg,ipg)
        dhxzdp = rlpxz_grid(irg,itg,ipg)
        dhyydp = rlpyy_grid(irg,itg,ipg)
        dhyzdp = rlpyz_grid(irg,itg,ipg)
        dhzzdp = rlpzz_grid(irg,itg,ipg)
        dhyxdp = dhxydp
        dhzxdp = dhxzdp
        dhzydp = dhyzdp

        trlie = gmxxu*dhxxdp + gmxyu*dhxydp + gmxzu*dhxzdp &
        &     + gmyxu*dhyxdp + gmyyu*dhyydp + gmyzu*dhyzdp &
        &     + gmzxu*dhzxdp + gmzyu*dhzydp + gmzzu*dhzzdp
!
        elpxx_grid(irg,itg,ipg) = dhxxdp - fa13*gmxxd*trlie
        elpxy_grid(irg,itg,ipg) = dhxydp - fa13*gmxyd*trlie
        elpxz_grid(irg,itg,ipg) = dhxzdp - fa13*gmxzd*trlie
        elpyy_grid(irg,itg,ipg) = dhyydp - fa13*gmyyd*trlie
        elpyz_grid(irg,itg,ipg) = dhyzdp - fa13*gmyzd*trlie
        elpzz_grid(irg,itg,ipg) = dhzzdp - fa13*gmzzd*trlie
!
      end do
    end do
  end do
!
end subroutine liegmab_gridpoint
