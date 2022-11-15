subroutine source_komar_mass_compact_qeos(souf)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrf, ntf, npf
  use def_matter
  use def_matter_parameter, only  :   radi, ber, ome
  use def_metric, only  :   psi, alph, bvxd, bvyd, bvzd
  use coordinate_grav_r, only       : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: souf(:,:,:)
  real(long), pointer :: alphf(:,:,:),  psif(:,:,:) 
  real(long), pointer :: bvxdf(:,:,:), bvydf(:,:,:), bvzdf(:,:,:)
  real(long)  ::    alphw, psiw, rhow, prew, hhw, utw, omew
  real(long)  ::   rjjx, rjjy, rjjz, rjjbeta, ene, rhoHw, esseS
  real(long)  ::   vphif(1:3)
  real(long)  ::   otermx, otermy, otermz, bvxufw, bvyufw, bvzufw, dummy
  integer     ::   irf,itf,ipf
  real(long)  ::   zfac, small = 1.0d-15
!
  call alloc_array3d(psif, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvxdf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvydf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvzdf, 0, nrf, 0, ntf, 0, npf)
!
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
  call interpo_gr2fl(bvxd, bvxdf)
  call interpo_gr2fl(bvyd, bvydf)
  call interpo_gr2fl(bvzd, bvzdf)
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        rhow = rhof(irf,itf,ipf)
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        bvxufw = bvxdf(irf,itf,ipf)
        bvyufw = bvydf(irf,itf,ipf)
        bvzufw = bvzdf(irf,itf,ipf)
        call quark_rho2phenedpdrho(rhow, prew, hhw, ene, dummy)
        utw  = utf(irf,itf,ipf)
        omew = omef(irf,itf,ipf)
!
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
!
        otermx = bvxufw + omew*vphif(1)
        otermy = bvyufw + omew*vphif(2)
        otermz = bvzufw + omew*vphif(3)
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
  deallocate(alphf)
  deallocate(psif)
  deallocate(bvxdf)
  deallocate(bvydf)
  deallocate(bvzdf)
end subroutine source_komar_mass_compact_qeos
