subroutine calc_h_udf_axisym
  use phys_constant, only  : long
  use grid_parameter
  use coordinate_grav_r, only : rg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_matter, only : rs, emd, hhf, utf , uxf , uyf , uzf, &
  &                                    utdf, uxdf, uydf, uzdf
  use def_matter_parameter, only : ome, ber
  use def_metric_on_SFC_CF, only : alphf, psif, bvxdf, bvydf, bvzdf, &
  &                                             bvxuf, bvyuf, bvzuf
  use def_metric_on_SFC_WL, only : hxxdf, hxydf, hxzdf, hyydf, hyzdf, hzzdf
  use def_emfield, only : vayd
  use def_vector_phi, only : vec_phif 
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use integrability_fnc_MHD
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_flgrad_4th_gridpoint
  implicit none
  real(long) :: hh, utd, ut, ux, uy, uz, vx, vy, vz, pre, rho, ene, qq
  real(long) :: bvdfc(3), ovufc(3), ovdfc(3)
  real(long) :: gmxxdf, gmxydf, gmxzdf, gmyxdf, gmyydf, gmyzdf, &
  &             gmzxdf, gmzydf, gmzzdf, ovovf,  alphff, psiff
  integer    :: irf, itf, ipf
!
  ipf = 0
  do itf = 0, ntf
    do irf = 0, nrf
!
      psiff  = psif(irf,itf,ipf)
      alphff = alphf(irf,itf,ipf)
      gmxxdf = 1.0d0 + hxxdf(irf,itf,ipf)
      gmxydf =         hxydf(irf,itf,ipf)
      gmxzdf =         hxzdf(irf,itf,ipf)
      gmyydf = 1.0d0 + hyydf(irf,itf,ipf)
      gmyzdf =         hyzdf(irf,itf,ipf)
      gmzzdf = 1.0d0 + hzzdf(irf,itf,ipf)
      gmyxdf = gmxydf
      gmzxdf = gmxzdf
      gmzydf = gmyzdf
!
      ut = utf(irf,itf,ipf)
      ux = uxf(irf,itf,ipf)
      uy = uyf(irf,itf,ipf)
      uz = uzf(irf,itf,ipf)
      vx = ux/ut
      vy = uy/ut
      vz = uz/ut
      bvdfc(1) = bvxdf(irf,itf,ipf)
      bvdfc(2) = bvydf(irf,itf,ipf)
      bvdfc(3) = bvzdf(irf,itf,ipf)
      ovufc(1) = bvxuf(irf,itf,ipf) + vx
      ovufc(2) = bvyuf(irf,itf,ipf) + vy
      ovufc(3) = bvzuf(irf,itf,ipf) + vz
      ovdfc(1) = bvxdf(irf,itf,ipf) + gmxxdf*vx + gmxydf*vy + gmxzdf*vz
      ovdfc(2) = bvydf(irf,itf,ipf) + gmyxdf*vx + gmyydf*vy + gmyzdf*vz
      ovdfc(3) = bvzdf(irf,itf,ipf) + gmzxdf*vx + gmzydf*vy + gmzzdf*vz
!
      utd = ut*(-alphff**2 + psiff**4*(bvdfc(1)*ovufc(1) &
      &          + bvdfc(2)*ovufc(2) + bvdfc(3)*ovufc(3)))
!
      qq = emd(irf,itf,ipf)
      call peos_q2hprho(qq, hh, pre, rho, ene)
      hhf( irf,itf,ipf) = hh
      utdf(irf,itf,ipf) = utd
      uxdf(irf,itf,ipf) = ut*psiff**4*ovdfc(1)
      uydf(irf,itf,ipf) = ut*psiff**4*ovdfc(2)
      uzdf(irf,itf,ipf) = ut*psiff**4*ovdfc(3)
!
    end do
  end do
!
! Copy to phi /= 0 planes.
!
  do ipf = 1, npf
    do itf = 0, ntf
      do irf = 0, nrf
        hhf(irf,itf,ipf)  = hhf(irf,itf,0)
        utdf(irf,itf,ipf) = utdf(irf,itf,0)
        uxdf(irf,itf,ipf) = cosphig(ipf)*uxdf(irf,itf,0) &
        &                 - sinphig(ipf)*uydf(irf,itf,0)
        uydf(irf,itf,ipf) = sinphig(ipf)*uxdf(irf,itf,0) &
        &                 + cosphig(ipf)*uydf(irf,itf,0)
        uzdf(irf,itf,ipf) = uzdf(irf,itf,0)
      end do
    end do
  end do
!
      itf = ntgeq; ipf = 0
!!      itf = 0; ipf = 0
      open(15,file='test_vec_4v',status='unknown')
        do irf = 0, nrf
          write(15,'(1p,9e20.12)')  rg(irf), utdf(irf,itf,ipf) &
              &                            , uxdf(irf,itf,ipf) &
              &                            , uydf(irf,itf,ipf) &
              &                            , uzdf(irf,itf,ipf)
        end do
      close(15)
!
end subroutine calc_h_udf_axisym
