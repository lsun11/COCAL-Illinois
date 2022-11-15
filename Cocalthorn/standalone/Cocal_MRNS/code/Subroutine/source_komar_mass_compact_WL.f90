subroutine source_komar_mass_compact_WL(souf)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only       : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_matter
  use def_matter_parameter, only  :   radi, ber, ome
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd, bvxu, bvyu, bvzu
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_metric_on_SFC_CF
  use def_metric_on_SFC_WL
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: souf(:,:,:)
  real(long) :: emdw, alphw, psiw, rhow, prew, hhw, utw, rhoHw, esseS
  real(long) :: rjjx, rjjy, rjjz, rjjbeta, ene
  real(long) :: vphif(1:3)
  real(long) :: otermx, otermy, otermz, omew
  real(long) :: bvxdfw, bvydfw, bvzdfw, bvxufw, bvyufw, bvzufw
  integer    :: irf, itf, ipf
  real(long) :: zfac, small = 1.0d-15
  real(long) :: hhxxdf, hhxydf, hhxzdf, hhyxdf, hhyydf, hhyzdf, &
  &             hhzxdf, hhzydf, hhzzdf
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        emdw = emd(irf,itf,ipf)
        if (emdw <= small) emdw = small
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        bvxufw = bvxuf(irf,itf,ipf)
        bvyufw = bvyuf(irf,itf,ipf)
        bvzufw = bvzuf(irf,itf,ipf)
        bvxdfw = bvxdf(irf,itf,ipf)
        bvydfw = bvydf(irf,itf,ipf)
        bvzdfw = bvzdf(irf,itf,ipf)
        call peos_q2hprho(emdw, hhw, prew, rhow, ene)
        utw = utf(irf,itf,ipf)
        omew = omef(irf,itf,ipf)
!
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
!
        hhxxdf = hxxdf(irf,itf,ipf)
        hhxydf = hxydf(irf,itf,ipf)
        hhxzdf = hxzdf(irf,itf,ipf)
        hhyydf = hyydf(irf,itf,ipf)
        hhyzdf = hyzdf(irf,itf,ipf)
        hhzzdf = hzzdf(irf,itf,ipf)
        hhyxdf = hhxydf
        hhzxdf = hhxzdf
        hhzydf = hhyzdf
!        
        otermx = bvxdfw + omew*vphif(1) &
        &      +   hhxxdf*omew*vphif(1) &
        &      +   hhxydf*omew*vphif(2) &
        &      +   hhxzdf*omew*vphif(3)
        otermy = bvydfw + omew*vphif(2) &
        &      +   hhyxdf*omew*vphif(1) &
        &      +   hhyydf*omew*vphif(2) &
        &      +   hhyzdf*omew*vphif(3)
        otermz = bvzdfw + omew*vphif(3) &
        &      +   hhzxdf*omew*vphif(1) &
        &      +   hhzydf*omew*vphif(2) &
        &      +   hhzzdf*omew*vphif(3)
!
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
end subroutine source_komar_mass_compact_WL
