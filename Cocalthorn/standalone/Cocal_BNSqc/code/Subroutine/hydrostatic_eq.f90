subroutine hydrostatic_eq(emd)
  use phys_constant, only  :   long
  use grid_parameter
  use def_matter, only : rs
  use def_matter_parameter, only : pinx, ome, ber
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_vector_phi, only : vec_phif
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: emd(:,:,:)
  real(long), pointer :: alphf(:,:,:), psif(:,:,:) 
  real(long), pointer :: bvxdf(:,:,:), bvydf(:,:,:), bvzdf(:,:,:)
  real(long) :: vphif(3)
  real(long) :: ovufc(3)
  real(long) :: ovufc2
  real(long) :: hh, ut
  integer    :: irf, itf, ipf
  call alloc_array3d(psif, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvxdf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvydf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(bvzdf, 0, nrf, 0, ntf, 0, npf)
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
  call interpo_gr2fl(bvxd, bvxdf)
  call interpo_gr2fl(bvyd, bvydf)
  call interpo_gr2fl(bvzd, bvzdf)
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        ovufc(1) = bvxdf(irf,itf,ipf) + ome*vphif(1)
        ovufc(2) = bvydf(irf,itf,ipf) + ome*vphif(2)
        ovufc(3) = bvzdf(irf,itf,ipf) + ome*vphif(3)
        ovufc2 = ovufc(1)**2 + ovufc(2)**2 + ovufc(3)**2
        ut = 1.0d0/sqrt(alphf(irf,itf,ipf)**2 & 
      &                - psif(irf,itf,ipf)**4*ovufc2)
        hh = ber*ut
        emd(irf,itf,ipf) = 1.0d0/(pinx+1.0d0)*(hh-1.0d0)
      end do
    end do
  end do
!
  deallocate(alphf)
  deallocate(psif)
  deallocate(bvxdf)
  deallocate(bvydf)
  deallocate(bvzdf)
end subroutine hydrostatic_eq
