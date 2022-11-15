subroutine hydrostatic_eq_qeos(rho)
  use phys_constant, only  :   long
  use grid_parameter
  use def_matter, only : utf, omef, jomef, jomef_int
  use def_matter_parameter, only : ome, ber, rhoc_qs, rhos_qs
  use def_metric_on_SFC_CF, only : alphf, psif, bvxdf, bvydf, bvzdf
  use coordinate_grav_r, only : rg
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: rho(:,:,:)
  real(long) :: vphif(3)
  real(long) :: ovufc(3), ovufc2
  real(long) :: omefc, jomef_intfc
  real(long) :: hh, ut, pre, ene, qq
  real(8)    :: rhomin, rhomax, rhotmp
  real(8), external :: quark_rho2h, quark_rho2h_dot
  integer    :: irf, itf, ipf
!
  rhomin = 0.5d0*rhos_qs
  rhomax = 1.1d0*rhoc_qs
  do ipf = 0, npf
    do itf = 0, ntf
      rhotmp = rhoc_qs
      do irf = 0, nrf
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        omefc    = omef(irf,itf,ipf)
        ovufc(1) = bvxdf(irf,itf,ipf) + omefc*vphif(1)
        ovufc(2) = bvydf(irf,itf,ipf) + omefc*vphif(2)
        ovufc(3) = bvzdf(irf,itf,ipf) + omefc*vphif(3)
        ovufc2 = ovufc(1)**2 + ovufc(2)**2 + ovufc(3)**2
        ut = 1.0d0/sqrt(alphf(irf,itf,ipf)**2 & 
      &                - psif(irf,itf,ipf)**4*ovufc2)
        jomef_intfc = jomef_int(irf,itf,ipf)
        hh = ber*ut*dexp(-jomef_intfc)
        call quark_h2rho(quark_rho2h,quark_rho2h_dot,rhomin,rhomax,hh,rhotmp) 
        rho(irf,itf,ipf) = rhotmp
        utf(irf,itf,ipf) = ut
      end do
    end do
  end do
!
end subroutine hydrostatic_eq_qeos
