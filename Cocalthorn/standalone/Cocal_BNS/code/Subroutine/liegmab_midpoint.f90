subroutine liegmab_midpoint
  use grid_parameter, only : nrg, ntg, npg
  use def_Lie_derivatives, only : elpxx, elpxy, elpxz, elpyy, elpyz, elpzz, &
  &                               rlpxx, rlpxy, rlpxz, rlpyy, rlpyz, rlpzz
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_dvphi, only : dphiu
  use interface_interpo_linear_type0
  use interface_grdphi_midpoint_type0
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
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(gmxxd,hxxd,irg,itg,ipg)
        call interpo_linear_type0(gmxyd,hxyd,irg,itg,ipg)
        call interpo_linear_type0(gmxzd,hxzd,irg,itg,ipg)
        call interpo_linear_type0(gmyyd,hyyd,irg,itg,ipg)
        call interpo_linear_type0(gmyzd,hyzd,irg,itg,ipg)
        call interpo_linear_type0(gmzzd,hzzd,irg,itg,ipg)
        gmxxd = gmxxd + 1.0d0
        gmyyd = gmyyd + 1.0d0
        gmzzd = gmzzd + 1.0d0
        gmyxd = gmxyd
        gmzxd = gmxzd
        gmzyd = gmyzd
        call interpo_linear_type0(gmxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(gmxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(gmxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(gmyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(gmyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(gmzzu,hzzu,irg,itg,ipg)
        gmxxu = gmxxu + 1.0d0
        gmyyu = gmyyu + 1.0d0
        gmzzu = gmzzu + 1.0d0
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu
!
        call grdphi_midpoint_type0(hxxd,dhxxdp,irg,itg,ipg)
        call grdphi_midpoint_type0(hxyd,dhxydp,irg,itg,ipg)
        call grdphi_midpoint_type0(hxzd,dhxzdp,irg,itg,ipg)
        call grdphi_midpoint_type0(hyyd,dhyydp,irg,itg,ipg)
        call grdphi_midpoint_type0(hyzd,dhyzdp,irg,itg,ipg)
        call grdphi_midpoint_type0(hzzd,dhzzdp,irg,itg,ipg)
!
        rlpxx(irg,itg,ipg) = dhxxdp + gmxyd*dphiu(2,1) + gmyxd*dphiu(2,1)
        rlpxy(irg,itg,ipg) = dhxydp + gmxxd*dphiu(1,2) + gmyyd*dphiu(2,1)
        rlpxz(irg,itg,ipg) = dhxzdp + gmyzd*dphiu(2,1)
        rlpyy(irg,itg,ipg) = dhyydp + gmyxd*dphiu(1,2) + gmxyd*dphiu(1,2)
        rlpyz(irg,itg,ipg) = dhyzdp + gmxzd*dphiu(1,2)
        rlpzz(irg,itg,ipg) = dhzzdp
!
        dhxxdp = rlpxx(irg,itg,ipg)
        dhxydp = rlpxy(irg,itg,ipg)
        dhxzdp = rlpxz(irg,itg,ipg)
        dhyydp = rlpyy(irg,itg,ipg)
        dhyzdp = rlpyz(irg,itg,ipg)
        dhzzdp = rlpzz(irg,itg,ipg)
        dhyxdp = dhxydp
        dhzxdp = dhxzdp
        dhzydp = dhyzdp
!
        trlie = gmxxu*dhxxdp + gmxyu*dhxydp + gmxzu*dhxzdp &
        &     + gmyxu*dhyxdp + gmyyu*dhyydp + gmyzu*dhyzdp &
        &     + gmzxu*dhzxdp + gmzyu*dhzydp + gmzzu*dhzzdp
!
        elpxx(irg,itg,ipg) = dhxxdp - fa13*gmxxd*trlie
        elpxy(irg,itg,ipg) = dhxydp - fa13*gmxyd*trlie
        elpxz(irg,itg,ipg) = dhxzdp - fa13*gmxzd*trlie
        elpyy(irg,itg,ipg) = dhyydp - fa13*gmyyd*trlie
        elpyz(irg,itg,ipg) = dhyzdp - fa13*gmyzd*trlie
        elpzz(irg,itg,ipg) = dhzzdp - fa13*gmzzd*trlie
!
      end do
    end do
  end do
!
end subroutine liegmab_midpoint
