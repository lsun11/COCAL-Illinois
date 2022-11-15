subroutine calc_4velocity_ut_spin
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF
  use def_matter
  use def_matter_parameter, only : ome, ber, confpow
  use def_vector_phi, only : hvec_phif, vec_phif
  use def_velocity_rot
  use def_velocity_potential
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use make_array_3d
  use interface_interpo_fl2gr
  implicit none
  real(long) :: psifc, alpfc
  real(long) :: emdfc, utfc, hhfc, prefc, rhofc, enefc
  real(long) :: vphif(3), ovdfc(3), ovdfc2, wx, wy, wz, psifc4, psifcp
  real(long) :: dxvep, dyvep, dzvep, lam, alpfc2, hut, dvep2, wdvep, w2, wterm, uih2
  integer :: ir, it, ip

!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        psifc = psif(ir,it,ip)
        alpfc= alphf(ir,it,ip)
        emdfc = emd(ir,it,ip)
!
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, enefc)
        vphif(1) = vec_phif(ir,it,ip,1)
        vphif(2) = vec_phif(ir,it,ip,2)
        vphif(3) = vec_phif(ir,it,ip,3)
        ovdfc(1) = bvxdf(ir,it,ip) + ome*vphif(1)
        ovdfc(2) = bvydf(ir,it,ip) + ome*vphif(2)
        ovdfc(3) = bvzdf(ir,it,ip) + ome*vphif(3)
        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,ir,it,ip)
        wx       = wxspf(ir,it,ip)
        wy       = wyspf(ir,it,ip)
        wz       = wzspf(ir,it,ip)

        psifc4   = psifc**4
        psifcp   = psifc**confpow
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

        utf(ir,it,ip) = hut/hhfc

      end do
    end do
  end do

  call interpo_fl2gr(utf,utg)
!
end subroutine calc_4velocity_ut_spin
