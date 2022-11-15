subroutine source_ang_mom(souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, rs, utf
  use def_matter_parameter
  use def_metric, only  :   tfkijkij, psi, alph, &
                            bvxd, bvyd, bvzd
  use coordinate_grav_r, only  :    rg
  use trigonometry_grav_theta, only  :   sinthg, costhg
  use trigonometry_grav_phi, only  :   sinphig, cosphig
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer ::   souf(:,:,:)
  integer     ::   ir,it,ip
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw
  real(long)  ::   rjjx, rjjy, rjjz, rjjphi
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   vphif(1:3)
  real(long)  ::   otermx, otermy, otermz, bvxufw, bvyufw, bvzufw
  real(long), pointer :: alphf(:,:,:),  psif(:,:,:) 
  real(long), pointer :: bvxdf(:,:,:), bvydf(:,:,:), bvzdf(:,:,:)
!
  call alloc_array3d(psif,  0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvxdf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvydf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvzdf, 0, nrf, 0, ntf, 0, npf)
!
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi , psif )
  call interpo_gr2fl(bvxd, bvxdf)
  call interpo_gr2fl(bvyd, bvydf)
  call interpo_gr2fl(bvzd, bvzdf)
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
        rhow = emdw**pinx
        prew = rhow*emdw
        hhw  = 1.0d0 + (pinx+1.0d0)*emdw
!
!        utw = hhw/ber
        utw = utf(ir,it,ip)
!
        vphif(1) = vec_phif(ir,it,ip,1)
        vphif(2) = vec_phif(ir,it,ip,2)
        vphif(3) = vec_phif(ir,it,ip,3)
!
        otermx = bvxufw + ome*vphif(1)
        otermy = bvyufw + ome*vphif(2)
        otermz = bvzufw + ome*vphif(3)
!
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
  deallocate(alphf)
  deallocate(psif)
  deallocate(bvxdf)
  deallocate(bvydf)
  deallocate(bvzdf)
end subroutine source_ang_mom
