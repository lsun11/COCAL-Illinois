subroutine source_ang_mom_peos_spin(souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, rs, vep, vepxf, vepyf, vepzf
  use def_metric_on_SFC_CF
  use def_velocity_rot
  use def_matter_parameter
  use def_metric, only  :   tfkijkij, psi, alph, &
                            bvxd, bvyd, bvzd
  use coordinate_grav_r, only  :    rg
  use trigonometry_grav_theta, only  :   sinthg, costhg
  use trigonometry_grav_phi, only  :   sinphig, cosphig
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer ::   souf(:,:,:)
  integer     ::   ir,it,ip
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, ene
  real(long)  ::   rjjx, rjjy, rjjz, rjjphi
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   vphif(1:3)
  real(long)  ::   otermx, otermy, otermz, bvxufw, bvyufw, bvzufw
  real(long)  ::   dxvep, dyvep, dzvep, lam, alphw2
  real(long)  ::   wx,wy,wz,psiw4,psiwp,dvep2,wdvep,w2,wterm,uih2,hut
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        emdw = emd(ir,it,ip)
        if (emdw <= small) emdw = small
        psiw = psif(ir,it,ip)
        alphw = alphf(ir,it,ip)
        bvxufw = bvxdf(ir,it,ip)
        bvyufw = bvydf(ir,it,ip)
        bvzufw = bvzdf(ir,it,ip)
        call peos_q2hprho(emdw, hhw, prew, rhow, ene)
!
        vphif(1) = vec_phif(ir,it,ip,1)
        vphif(2) = vec_phif(ir,it,ip,2)
        vphif(3) = vec_phif(ir,it,ip,3)
        otermx = bvxufw + ome*vphif(1)
        otermy = bvyufw + ome*vphif(2)
        otermz = bvzufw + ome*vphif(3)
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
        lam      = ber + otermx*dxvep + otermy*dyvep + otermz*dzvep
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

        zfac = 1.0d0
        if (emdw <= small) zfac = 0.0d0
        rjjx = hhw*rhow*alphw*utw**2*psiw**4*otermx
        rjjy = hhw*rhow*alphw*utw**2*psiw**4*otermy
        rjjz = hhw*rhow*alphw*utw**2*psiw**4*otermz
!
        rjjphi = rjjx*vphif(1) + rjjy*vphif(2) + rjjz*vphif(3)
        souf(ir,it,ip) = rjjphi*psiw**6*zfac
!
      end do
    end do
  end do
!
end subroutine source_ang_mom_peos_spin
