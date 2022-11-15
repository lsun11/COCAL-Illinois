subroutine calc_utg_WL
  use phys_constant, only  : long
  use grid_parameter
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd, bvxu, bvyu, bvzu
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_matter, only : utg, omeg, emdg
  use def_matter_parameter, only : ome, ber, emdc
  use coordinate_grav_r, only : rg
  use def_vector_phi, only : vec_phig
  use make_array_3d
  implicit none
  real(long) :: vphig(3)
  real(long) :: ovugc(3), ovdgc(3)
  real(long) :: omegc, jomeg_intgc
  real(long) :: hh, ut, pre, rho, ene, qq, zfac
  real(long) :: gmxxd, gmxyd, gmxzd, gmyxd, gmyyd, gmyzd, &
  &             gmzxd, gmzyd, gmzzd, ovovg
  integer    :: irg, itg, ipg
!
  zfac = 5.0d0*emdg(0,0,0)

  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        qq = emdg(irg,itg,ipg)

        if (qq > zfac) then
          vphig(1) = vec_phig(irg,itg,ipg,1)
          vphig(2) = vec_phig(irg,itg,ipg,2)
          vphig(3) = vec_phig(irg,itg,ipg,3)

          omegc    = omeg(irg,itg,ipg)

          gmxxd = 1.0d0 + hxxd(irg,itg,ipg)
          gmxyd =         hxyd(irg,itg,ipg)
          gmxzd =         hxzd(irg,itg,ipg)
          gmyyd = 1.0d0 + hyyd(irg,itg,ipg)
          gmyzd =         hyzd(irg,itg,ipg)
          gmzzd = 1.0d0 + hzzd(irg,itg,ipg)
          gmyxd = gmxyd
          gmzxd = gmxzd
          gmzyd = gmyzd

          ovugc(1) = bvxu(irg,itg,ipg) + omegc*vphig(1)
          ovugc(2) = bvyu(irg,itg,ipg) + omegc*vphig(2)
          ovugc(3) = bvzu(irg,itg,ipg) + omegc*vphig(3)

          ovdgc(1) = bvxd(irg,itg,ipg) &
          &        + gmxxd*omegc*vphig(1) + gmxyd*omegc*vphig(2) &
          &        + gmxzd*omegc*vphig(3)
          ovdgc(2) = bvyd(irg,itg,ipg) &
          &        + gmyxd*omegc*vphig(1) + gmyyd*omegc*vphig(2) &
          &        + gmyzd*omegc*vphig(3)
          ovdgc(3) = bvzd(irg,itg,ipg) &
          &        + gmzxd*omegc*vphig(1) + gmzyd*omegc*vphig(2) &
          &        + gmzzd*omegc*vphig(3)

          ovovg = ovdgc(1)*ovugc(1) + ovdgc(2)*ovugc(2) + ovdgc(3)*ovugc(3)

          ut = 1.0d0/sqrt(alph(irg,itg,ipg)**2 & 
          &                - psi(irg,itg,ipg)**4*ovovg)

          utg(irg,itg,ipg) = ut
        end if
      end do
    end do
  end do
!
end subroutine calc_utg_WL
