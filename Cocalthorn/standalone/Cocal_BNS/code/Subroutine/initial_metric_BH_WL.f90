subroutine initial_metric_BH_WL
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_metric
  use def_metric_hij
  use def_matter, only : omeg
  use def_metric_excurve_grid, only : trk_grid
  use def_vector_x, only : vec_xg, hvec_xg
  implicit none
  real(long) :: psi_kerr, alph_kerr, bvu_kerr(3), bvd_kerr(3), &
  &             hijd_kerr(3,3), hiju_kerr(3,3), x, y, z, traceK
  integer :: irg, itg, ipg
!
!!!  psi(0:nrg,0:ntg,0:npg)  = 1.0d0
!!!  alph(0:nrg,0:ntg,0:npg) = 1.0d0
!!!  alps(0:nrg,0:ntg,0:npg) = 1.0d0
!!!  bvxu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  bvyu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  bvzu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  bvxd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  bvyd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  bvzd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hxxd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hxyd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hxzd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hyyd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hyzd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hzzd(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hxxu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hxyu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hxzu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hyyu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hyzu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!  hzzu(0:nrg,0:ntg,0:npg) = 0.0d0
!!!
  call calc_vector_x_grav(1)
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        x = vec_xg(irg,itg,ipg,1)
        y = vec_xg(irg,itg,ipg,2)
        z = vec_xg(irg,itg,ipg,3)
        call kerr_schild_metric_3plus1(x,y,z,psi_kerr,alph_kerr, &
        &                              bvu_kerr,bvd_kerr,hijd_kerr,hiju_kerr)
        psi(irg,itg,ipg)  = psi_kerr
        alph(irg,itg,ipg) = alph_kerr
        alps(irg,itg,ipg) = psi_kerr*alph_kerr
        bvxu(irg,itg,ipg) = bvu_kerr(1)
        bvyu(irg,itg,ipg) = bvu_kerr(2)
        bvzu(irg,itg,ipg) = bvu_kerr(3)
        bvxd(irg,itg,ipg) = bvd_kerr(1)
        bvyd(irg,itg,ipg) = bvd_kerr(2)
        bvzd(irg,itg,ipg) = bvd_kerr(3)
        hxxd(irg,itg,ipg) = hijd_kerr(1,1)
        hxyd(irg,itg,ipg) = hijd_kerr(1,2)
        hxzd(irg,itg,ipg) = hijd_kerr(1,3)
        hyyd(irg,itg,ipg) = hijd_kerr(2,2)
        hyzd(irg,itg,ipg) = hijd_kerr(2,3)
        hzzd(irg,itg,ipg) = hijd_kerr(3,3)
        hxxu(irg,itg,ipg) = hiju_kerr(1,1)
        hxyu(irg,itg,ipg) = hiju_kerr(1,2)
        hxzu(irg,itg,ipg) = hiju_kerr(1,3)
        hyyu(irg,itg,ipg) = hiju_kerr(2,2)
        hyzu(irg,itg,ipg) = hiju_kerr(2,3)
        hzzu(irg,itg,ipg) = hiju_kerr(3,3)
        call kerr_schild_traceK(x,y,z,traceK)
        trk_grid(irg,itg,ipg) = traceK
        if (irg.ne.0.and.itg.ne.0.and.ipg.ne.0) then
          x = hvec_xg(irg,itg,ipg,1)
          y = hvec_xg(irg,itg,ipg,2)
          z = hvec_xg(irg,itg,ipg,3)
          call kerr_schild_traceK(x,y,z,traceK)
          trk(irg,itg,ipg) = traceK
        end if
      end do
    end do
  end do
  omeg(0:nrg,0:ntg,0:npg) = 0.0d0  ! For rotating shift
!
      itg = ntg/2-2; ipg = npg/2-2
      open(16,file='test_vec0',status='unknown')
        do irg = 0, nrg
          write(16,'(1p,13e20.12)') rg(irg),  psi(irg,itg,ipg)  &
              &                            , alph(irg,itg,ipg) &
              &                            , bvxd(irg,itg,ipg) &
              &                            , bvyd(irg,itg,ipg) &
              &                            , bvzd(irg,itg,ipg) &
              &                            , hxxd(irg,itg,ipg) &
              &                            , hxyd(irg,itg,ipg) &
              &                            , hxzd(irg,itg,ipg) &
              &                            , hyyd(irg,itg,ipg) &
              &                            , hyzd(irg,itg,ipg) &
              &                            , hzzd(irg,itg,ipg) &
              &                            , trk_grid(irg,itg,ipg)
        end do
      close(16)
!
end subroutine initial_metric_BH_WL
