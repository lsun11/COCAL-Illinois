subroutine source_rest_mass_peos_spin(souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, vep, vepxf, vepyf, vepzf
  use def_metric_on_SFC_CF
  use def_velocity_rot
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_matter_parameter, only  :   radi, ber, ome, confpow
  use def_metric, only  :   tfkijkij,  psi, alph
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_interpo_gr2fl
  implicit none
  real(long),pointer ::   souf(:,:,:)
  real(long)  ::   psiwm6
  real(long)  ::   emdw, alphw, psiw, rhow, hhw, utw, rhoHw, esseS
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   otermx, otermy, otermz, prew, ene
  real(long)  ::   dxvep, dyvep, dzvep, lam, alphw2, ovdfc(3)
  integer     ::   ir,it,ip
  real(long)  ::   wx,wy,wz,psiw4,psiwp,dvep2,wdvep,w2,wterm,uih2,hut
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        emdw = emd(ir,it,ip)
        if (emdw <= small) emdw = small
        psiw = psif(ir,it,ip)
        alphw = alphf(ir,it,ip)
        call peos_q2hprho(emdw, hhw, prew, rhow, ene)

        ovdfc(1) = bvxdf(ir,it,ip) + ome*vec_phif(ir,it,ip,1)
        ovdfc(2) = bvydf(ir,it,ip) + ome*vec_phif(ir,it,ip,2)
        ovdfc(3) = bvzdf(ir,it,ip) + ome*vec_phif(ir,it,ip,3)
        dxvep = vepxf(ir,it,ip) 
        dyvep = vepyf(ir,it,ip)
        dzvep = vepzf(ir,it,ip)
!
        wx       = wxspf(ir,it,ip)
        wy       = wyspf(ir,it,ip)
        wz       = wzspf(ir,it,ip)
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

        souf(ir,it,ip) = rhow*alphw*utw*psiw**6
      end do
    end do
  end do
!
end subroutine source_rest_mass_peos_spin
