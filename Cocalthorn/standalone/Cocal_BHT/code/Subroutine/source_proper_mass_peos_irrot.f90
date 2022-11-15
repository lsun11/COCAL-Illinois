subroutine source_proper_mass_peos_irrot(souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, vep, vepxf, vepyf, vepzf
  use def_metric_on_SFC_CF
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_matter_parameter, only  :   radi, ber, pinx, ome
  use def_metric, only  :   tfkijkij,  psi, alph
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_interpo_gr2fl
  implicit none
  real(long),pointer ::   souf(:,:,:)
  real(long)  ::   psiwm6
  real(long)  ::   emdw, alphw, psiw, epsilonw, hhw, utw, rhoHw, esseS
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   otermx, otermy, otermz, prew, rhow
  real(long)  ::   dxvep, dyvep, dzvep, lam, alphw2, ovdfc(3)
  integer     ::   ir,it,ip
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        emdw = emd(ir,it,ip)
        if (emdw <= small) emdw = small
        psiw = psif(ir,it,ip)
        alphw = alphf(ir,it,ip)
        call peos_q2hprho(emdw, hhw, prew, rhow, epsilonw)

        ovdfc(1) = bvxdf(ir,it,ip) + ome*vec_phif(ir,it,ip,1)
        ovdfc(2) = bvydf(ir,it,ip) + ome*vec_phif(ir,it,ip,2)
        ovdfc(3) = bvzdf(ir,it,ip) + ome*vec_phif(ir,it,ip,3)
        dxvep = vepxf(ir,it,ip)
        dyvep = vepyf(ir,it,ip)
        dzvep = vepzf(ir,it,ip)
!        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,ir,it,ip)
        alphw2 = alphw**2
        lam    = ber + ovdfc(1)*dxvep + ovdfc(2)*dyvep + ovdfc(3)*dzvep
        utw = lam/alphw2/hhw
!
        souf(ir,it,ip) = epsilonw*alphw*utw*psiw**6
      end do
    end do
  end do
!
end subroutine source_proper_mass_peos_irrot
