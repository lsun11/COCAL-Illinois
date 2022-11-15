subroutine hydrostatic_eq_CF_peos_lecc_irrot(emd,utf,vepxf,vepyf,vepzf)
  use phys_constant, only  :   long
  use grid_parameter
  use def_matter, only : rs, omef, jomef, jomef_int, &
  &                      omeg, jomeg, jomeg_int, vep
  use def_velocity_potential
  use def_matter_parameter, only : ome, ber, ROT_LAW, velx
  use def_metric_on_SFC_CF
  use coordinate_grav_r, only : rg
  use def_vector_phi, only : vec_phif
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

        ovdfc(1) = bvxdf(irf,itf,ipf) + omefc*vphif(1) + velx
        ovdfc(2) = bvydf(irf,itf,ipf) + omefc*vphif(2)
        ovdfc(3) = bvzdf(irf,itf,ipf) + omefc*vphif(3)

        call flgrad_2nd_gridpoint(vep,dxvep,dyvep,dzvep,irf,itf,ipf)

        psifc    = psif(irf,itf,ipf)
        psifc4   = psifc**4
        alpfc    = alphf(irf,itf,ipf)
        alpfc2   = alpfc**2
        lam      = ber + ovdfc(1)*dxvep + ovdfc(2)*dyvep + ovdfc(3)*dzvep
     
        dvep2    = (dxvep**2 + dyvep**2 + dzvep**2)/psifc4
        uih2     = dvep2 

        hut = lam/alpfc2

        if ( (hut*hut*alpfc2 - uih2)<0.0d0 ) then
          write(6,*)  "hh imaginary....exiting"
          stop
        end if
        hh = sqrt(hut*hut*alpfc2 - uih2)

!        utf(irf,itf,ipf) = hut/hh

!        vepxf(irf,itf,ipf) = dxvep
!        vepyf(irf,itf,ipf) = dyvep
!        vepzf(irf,itf,ipf) = dzvep

        call peos_h2qprho(hh, qq, pre, rho, ene)
        emd(irf,itf,ipf) = qq
      end do
    end do
  end do
!
end subroutine hydrostatic_eq_CF_peos_lecc_irrot
