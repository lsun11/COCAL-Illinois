subroutine hydrostatic_eq_CF_peos_spin(emd,utf,vepxf,vepyf,vepzf)
  use phys_constant, only  :   long
  use grid_parameter
  use def_matter, only : rs, omef, jomef, jomef_int, &
  &                      omeg, jomeg, jomeg_int, vep
  use def_velocity_potential
  use def_velocity_rot
  use def_matter_parameter, only : ome, ber, ROT_LAW,confpow
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use def_metric_on_SFC_CF
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_vector_phi, only : vec_phif
  use def_vector_x, only : vec_xf, vec_xg
  use make_array_3d
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  implicit none
  real(long), pointer :: emd(:,:,:), utf(:,:,:), vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
  real(long) :: vphif(3)
  real(long) :: ovdfc(3), ovdfc2
  real(long) :: dxvep, dyvep, dzvep, lam, wx, wy, wz, wterm
  real(long) :: psifc, psifc4, alpfc, alpfc2, hut, psifcp
  real(long) :: dvep2, wdvep, w2, uih2
  real(long) :: alp_tmp, psi_tmp, bx_tmp, by_tmp
  real(long) :: omefc, jomefc, jomef_intfc, omegc, jomegc, jomeg_intgc
  real(long) :: hh, ut, pre, rho, ene, qq
  real(long) :: xx, yy, zz, Rcyl
  integer    :: irf, itf, ipf
!
  if (ROT_LAW.eq.'DR') then 
    do itf = 0, ntf
      do irf = 0, nrf
        ipf = 0
        xx = vec_xf(irf,itf,ipf,1)
        yy = vec_xf(irf,itf,ipf,2)
        zz = vec_xf(irf,itf,ipf,3)
        Rcyl = sqrt(xx**2 + yy**2)
        omefc = omef(irf,itf,ipf)
!
        alp_tmp = alphf(irf,itf,ipf)
        psi_tmp = psif(irf,itf,ipf)
        bx_tmp = bvxdf(irf,itf,ipf)
        by_tmp = bvydf(irf,itf,ipf)
!        
        call calc_omega_drot(Rcyl,alp_tmp,psi_tmp,bx_tmp,by_tmp&
             & ,omefc,jomefc,jomef_intfc)
        omef(irf,itf,0:npf) = omefc
        jomef(irf,itf,0:npf) = jomefc
        jomef_int(irf,itf,0:npf) = jomef_intfc
!
        xx = vec_xg(irf,itf,ipf,1)
        yy = vec_xg(irf,itf,ipf,2)
        zz = vec_xg(irf,itf,ipf,3)
        Rcyl = sqrt(xx**2 + yy**2)
        omegc = omeg(irf,itf,ipf)
        alp_tmp = alph(irf,itf,ipf)
        psi_tmp = psi(irf,itf,ipf)
        bx_tmp = bvxd(irf,itf,ipf)
        by_tmp = bvyd(irf,itf,ipf)
        call calc_omega_drot(Rcyl,alp_tmp,psi_tmp,bx_tmp,by_tmp&
             & ,omegc,jomegc,jomeg_intgc)
        omeg(irf,itf,0:npf) = omegc
        jomeg(irf,itf,0:npf) = jomegc
        jomeg_int(irf,itf,0:npf) = jomeg_intgc
      end do
    end do
  else 
    omef(0:nrf,0:ntf,0:npf) = ome
    jomef(0:nrf,0:ntf,0:npf) = 0.0d0
    jomef_int(0:nrf,0:ntf,0:npf) = 0.0d0
    omeg(0:nrf,0:ntf,0:npf) = ome
    jomeg(0:nrf,0:ntf,0:npf) = 0.0d0
    jomeg_int(0:nrf,0:ntf,0:npf) = 0.0d0
  end if
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        omefc    = omef(irf,itf,ipf)
        ovdfc(1) = bvxdf(irf,itf,ipf) + omefc*vphif(1)
        ovdfc(2) = bvydf(irf,itf,ipf) + omefc*vphif(2)
        ovdfc(3) = bvzdf(irf,itf,ipf) + omefc*vphif(3)
        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,irf,itf,ipf)
!        call flgrad_4th_gridpoint(vep,dxvep,dyvep,dzvep,irf,itf,ipf)
        wx       = wxspf(irf,itf,ipf)
        wy       = wyspf(irf,itf,ipf)
        wz       = wzspf(irf,itf,ipf)
        psifc    = psif(irf,itf,ipf)
        psifc4   = psifc**4
        psifcp   = psifc**confpow
        alpfc    = alphf(irf,itf,ipf)
        alpfc2   = alpfc**2
        lam      = ber + ovdfc(1)*dxvep + ovdfc(2)*dyvep + ovdfc(3)*dzvep
     
        dvep2    = (dxvep**2 + dyvep**2 + dzvep**2)/psifc4
        wdvep    = (wx*dxvep + wy*dyvep + wz*dzvep)*psifcp
        w2       = psifc4*(wx*wx + wy*wy + wz*wz)*psifcp**2.0d0

        wterm    = wdvep + w2
        uih2     = dvep2 + 2.0d0*wdvep + w2

        if ( (lam*lam + 4.0d0*alpfc2*wterm)<0.0d0 ) then
          write(6,*)  "hut imaginary....exiting"
          stop
        end if  
        hut = (lam + sqrt(lam*lam + 4.0d0*alpfc2*wterm))/(2.0d0*alpfc2)

        if ( (hut*hut*alpfc2 - uih2)<0.0d0 ) then
          write(6,*)  "hh imaginary....exiting"
          stop
        end if
        hh = sqrt(hut*hut*alpfc2 - uih2)

!        utf(irf,itf,ipf)  = hut/hh

        call peos_h2qprho(hh, qq, pre, rho, ene)
        emd(irf,itf,ipf) = qq
      end do
    end do
  end do
!
end subroutine hydrostatic_eq_CF_peos_spin
