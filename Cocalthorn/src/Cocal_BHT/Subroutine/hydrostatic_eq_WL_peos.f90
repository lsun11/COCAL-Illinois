subroutine hydrostatic_eq_WL_peos(emd)
  use phys_constant, only  : long
  use grid_parameter
  use def_matter, only : utf, omef, jomef, jomef_int
  use def_matter_parameter, only : ome, ber
  use def_metric_on_SFC_CF
  use def_metric_on_SFC_WL
  use coordinate_grav_r, only : rg
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: emd(:,:,:)
  real(long) :: vphif(3)
  real(long) :: ovufc(3), ovdfc(3)
  real(long) :: omefc, jomef_intfc
  real(long) :: hh, ut, pre, rho, ene, qq
  real(long) :: gmxxdf, gmxydf, gmxzdf, gmyxdf, gmyydf, gmyzdf, &
  &             gmzxdf, gmzydf, gmzzdf, ovovf
  integer    :: irf, itf, ipf
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        omefc    = omef(irf,itf,ipf)
        gmxxdf = 1.0d0 + hxxdf(irf,itf,ipf)
        gmxydf =         hxydf(irf,itf,ipf)
        gmxzdf =         hxzdf(irf,itf,ipf)
        gmyydf = 1.0d0 + hyydf(irf,itf,ipf)
        gmyzdf =         hyzdf(irf,itf,ipf)
        gmzzdf = 1.0d0 + hzzdf(irf,itf,ipf)
        gmyxdf = gmxydf
        gmzxdf = gmxzdf
        gmzydf = gmyzdf
        ovufc(1) = bvxuf(irf,itf,ipf) + omefc*vphif(1)
        ovufc(2) = bvyuf(irf,itf,ipf) + omefc*vphif(2)
        ovufc(3) = bvzuf(irf,itf,ipf) + omefc*vphif(3)
        ovdfc(1) = bvxdf(irf,itf,ipf) &
        &        + gmxxdf*omefc*vphif(1) + gmxydf*omefc*vphif(2) &
        &        + gmxzdf*omefc*vphif(3)
        ovdfc(2) = bvydf(irf,itf,ipf) &
        &        + gmyxdf*omefc*vphif(1) + gmyydf*omefc*vphif(2) &
        &        + gmyzdf*omefc*vphif(3)
        ovdfc(3) = bvzdf(irf,itf,ipf) &
        &        + gmzxdf*omefc*vphif(1) + gmzydf*omefc*vphif(2) &
        &        + gmzzdf*omefc*vphif(3)
        ovovf = ovdfc(1)*ovufc(1) + ovdfc(2)*ovufc(2) + ovdfc(3)*ovufc(3)
        ut = 1.0d0/sqrt(alphf(irf,itf,ipf)**2 & 
      &                - psif(irf,itf,ipf)**4*ovovf)
        jomef_intfc = jomef_int(irf,itf,ipf)
        hh = ber*ut*dexp(-jomef_intfc)
        call peos_h2qprho(hh, qq, pre, rho, ene)
        emd(irf,itf,ipf) = qq
        utf(irf,itf,ipf) = ut
      end do
    end do
  end do
!
end subroutine hydrostatic_eq_WL_peos
