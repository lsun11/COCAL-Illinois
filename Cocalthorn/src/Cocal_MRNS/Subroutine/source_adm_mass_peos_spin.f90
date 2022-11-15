subroutine source_adm_mass_peos_spin(soug,souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, vep, vepxf, vepyf, vepzf
  use def_metric_on_SFC_CF
  use def_velocity_rot
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_matter_parameter, only  :  radi, ber, ome, confpow
  use def_metric, only  :   tfkijkij, psi, alph
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_interpo_gr2fl
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: soug(:,:,:)
  real(long), pointer :: souf(:,:,:)
  integer     ::   irg,itg,ipg,irf,itf,ipf
  real(long)  ::   psiwm7
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, rhoHw
  real(long)  ::   epsilonw, alutw
  real(long)  ::   dxvep, dyvep, dzvep, lam, alphw2, ovdfc(3)
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   wx,wy,wz,psiw4,psiwp,dvep2,wdvep,w2,wterm,uih2,hut
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        soug(irg,itg,ipg) = 0.125d0*psiw**5*tfkijkij(irg,itg,ipg)
      end do
    end do
  end do
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        emdw = emd(irf,itf,ipf)
        if (emdw <= small) emdw = small
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        call peos_q2hprho(emdw, hhw, prew, rhow, epsilonw)

        ovdfc(1) = bvxdf(irf,itf,ipf) + ome*vec_phif(irf,itf,ipf,1)
        ovdfc(2) = bvydf(irf,itf,ipf) + ome*vec_phif(irf,itf,ipf,2)
        ovdfc(3) = bvzdf(irf,itf,ipf) + ome*vec_phif(irf,itf,ipf,3)
        dxvep = vepxf(irf,itf,ipf)
        dyvep = vepyf(irf,itf,ipf)
        dzvep = vepzf(irf,itf,ipf)
!
        wx       = wxspf(irf,itf,ipf)
        wy       = wyspf(irf,itf,ipf)
        wz       = wzspf(irf,itf,ipf)
        psiw4    = psiw**4
        psiwp    = psiw**confpow
        alphw2   = alphw**2
        lam      = ber + ovdfc(1)*dxvep + ovdfc(2)*dyvep + ovdfc(3)*dzvep
        dvep2    = (dxvep**2 + dyvep**2 + dzvep**2)/psiw4
        wdvep    = (wx*dxvep + wy*dyvep + wz*dzvep)*psiwp
        w2       = psiw4*(wx*wx + wy*wy + wz*wz)*psiwp**2.0d0
        wterm    = wdvep + w2
        uih2     = dvep2 + 2.0d0*wdvep + w2
        if ( (lam*lam + 4.0d0*alphw2*wterm)<0.0d0 ) then
          write(6,*)  "hut imaginary....exiting"
          stop
        end if
        hut = (lam + sqrt(lam*lam + 4.0d0*alphw2*wterm))/(2.0d0*alphw2)
        utw = hut/hhw
        alutw = alphw*utw
!
!        souf(irf,itf,ipf) = 2.0d0*pi*psiw**5*rhoHw
        souf(irf,itf,ipf) = 2.0d0*pi*psiw**5  &
        &  *(epsilonw*alutw**2 + prew*(alutw-1.0d0)*(alutw+1.0d0))
      end do
    end do
  end do
!
end subroutine source_adm_mass_peos_spin
