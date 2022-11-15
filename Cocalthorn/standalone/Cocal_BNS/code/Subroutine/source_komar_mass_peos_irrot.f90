subroutine source_komar_mass_peos_irrot(soug,souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, vep, vepxf, vepyf, vepzf
  use def_metric_on_SFC_CF
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_matter_parameter, only  :   radi, ber, ome
  use def_metric, only  :   tfkijkij,  psi, alph
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_interpo_gr2fl
  use interface_interpo_linear_type0
  implicit none
  real(long),pointer ::   soug(:,:,:), souf(:,:,:)
  real(long)  ::   psiwm6
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, rhoHw, esseS
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   otermx, otermy, otermz, ene
  real(long)  ::   dxvep, dyvep, dzvep, lam, alphw2, ovdfc(3)
  integer     ::   irg,itg,ipg,irf,itf,ipf
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        call interpo_linear_type0(alphw,alph,irg,itg,ipg)
        psiwm6 = psiw**6
        soug(irg,itg,ipg) = alphw*psiwm6*tfkijkij(irg,itg,ipg)
      end do
    end do
  end do
!
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        emdw = emd(irf,itf,ipf)
        if (emdw <= small) emdw = small
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        call peos_q2hprho(emdw, hhw, prew, rhow, ene)

        ovdfc(1) = bvxdf(irf,itf,ipf) + ome*vec_phif(irf,itf,ipf,1)
        ovdfc(2) = bvydf(irf,itf,ipf) + ome*vec_phif(irf,itf,ipf,2)
        ovdfc(3) = bvzdf(irf,itf,ipf) + ome*vec_phif(irf,itf,ipf,3)
        dxvep = vepxf(irf,itf,ipf)
        dyvep = vepyf(irf,itf,ipf)
        dzvep = vepzf(irf,itf,ipf)

!        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,irf,itf,ipf)
        alphw2 = alphw**2
        lam    = ber + ovdfc(1)*dxvep + ovdfc(2)*dyvep + ovdfc(3)*dzvep
        utw = lam/alphw2/hhw
!
        rhoHw = hhw*rhow*(alphw*utw)**2 - prew
        esseS = -hhw*rhow + 4.0d0*prew + rhoHw
! 
        souf(irf,itf,ipf) = 4.0d0*pi*alphw*psiw**6*(esseS+rhoHw)
      end do
    end do
  end do
!
end subroutine source_komar_mass_peos_irrot
