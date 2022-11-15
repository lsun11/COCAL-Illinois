subroutine source_komar_mass_compact_peos_irrot(souf)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_matter
  use def_matter_parameter, only  :   radi, ber, ome
  use def_metric, only  :   psi, alph, bvxd, bvyd, bvzd
  use coordinate_grav_r, only       : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: souf(:,:,:)
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, rhoHw, esseS
  real(long)  ::   rjjx, rjjy, rjjz, rjjbeta, ene
  real(long)  ::   vphif(1:3)
  real(long)  ::   otermx, otermy, otermz, bvxufw, bvyufw, bvzufw
  integer     ::   irf,itf,ipf
  real(long)  ::   dxvep, dyvep, dzvep, lam, alphw2
  real(long)  ::   zfac, small = 1.0d-15
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        emdw = emd(irf,itf,ipf)
        if (emdw <= small) emdw = small
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        bvxufw = bvxdf(irf,itf,ipf)
        bvyufw = bvydf(irf,itf,ipf)
        bvzufw = bvzdf(irf,itf,ipf)
        call peos_q2hprho(emdw, hhw, prew, rhow, ene)
!
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        otermx = bvxufw + ome*vphif(1)
        otermy = bvyufw + ome*vphif(2)
        otermz = bvzufw + ome*vphif(3)
        dxvep = vepxf(irf,itf,ipf)
        dyvep = vepyf(irf,itf,ipf)
        dzvep = vepzf(irf,itf,ipf)
!        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,irf,itf,ipf)
        alphw2 = alphw**2
        lam    = ber + otermx*dxvep + otermy*dyvep + otermz*dzvep
        utw = lam/alphw2/hhw

        rhoHw = hhw*rhow*(alphw*utw)**2 - prew
        esseS = -hhw*rhow + 4.0d0*prew + rhoHw
! 
        rjjx = hhw*rhow*alphw*utw**2*psiw**4*otermx
        rjjy = hhw*rhow*alphw*utw**2*psiw**4*otermy
        rjjz = hhw*rhow*alphw*utw**2*psiw**4*otermz
!
        rjjbeta = rjjx*bvxufw + rjjy*bvyufw + rjjz*bvzufw
!
        souf(irf,itf,ipf) = (alphw*(esseS+rhoHw) - 2.0d0*rjjbeta)*psiw**6
      end do
    end do
  end do
!
end subroutine source_komar_mass_compact_peos_irrot
