subroutine hydrostatic_eq_peos_lecc(emd)
  use phys_constant, only  :   long
  use grid_parameter
  use def_matter, only : utf, omef, jomef, jomef_int
  use def_matter_parameter, only : ome, ber, velx
  use def_metric_on_SFC_CF, only : alphf, psif, bvxdf, bvydf, bvzdf
  use coordinate_grav_r, only : rg
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: emd(:,:,:)
  real(long) :: vphif(3)
  real(long) :: ovufc(3), ovufc2
  real(long) :: omefc, jomef_intfc
  real(long) :: hh, ut, pre, rho, ene, qq
  integer    :: irf, itf, ipf
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        omefc    = omef(irf,itf,ipf)
        ovufc(1) = bvxdf(irf,itf,ipf) + omefc*vphif(1) + velx
        ovufc(2) = bvydf(irf,itf,ipf) + omefc*vphif(2)
        ovufc(3) = bvzdf(irf,itf,ipf) + omefc*vphif(3)
        ovufc2 = ovufc(1)**2 + ovufc(2)**2 + ovufc(3)**2
        ut = 1.0d0/sqrt(alphf(irf,itf,ipf)**2 & 
      &                - psif(irf,itf,ipf)**4*ovufc2)
        jomef_intfc = jomef_int(irf,itf,ipf)
        hh = ber*ut*dexp(-jomef_intfc)
        call peos_h2qprho(hh, qq, pre, rho, ene)
        emd(irf,itf,ipf) = qq
        utf(irf,itf,ipf) = ut
      end do
    end do
  end do
!
end subroutine hydrostatic_eq_peos_lecc
