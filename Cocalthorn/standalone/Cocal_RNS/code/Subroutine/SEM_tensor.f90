subroutine SEM_tensor
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF, only : psif, alphf, bvxdf, bvydf, bvzdf
  use def_metric_on_SFC_WL, only : hxxdf, hxydf, hxzdf, hyydf, hyzdf, hzzdf, &
  &                                hxxuf, hxyuf, hxzuf, hyyuf, hyzuf, hzzuf
  use def_matter,     only : emd , utf, uxf, uyf, uzf
  use def_SEM_tensor, only : rhoH, jmd, smijd, trsm
  implicit none
  integer :: irf, itf, ipf, ia, ib
  real(long) :: psifc, psifc4, psifc4inv, psifc8, &
  &             alpfc, bvxdfc, bvydfc, bvzdfc
  real(long) :: emdfc, hhfc, prefc, rhofc, enefc
  real(long) :: utfc, uxfc, uyfc, uzfc, ut, ux, uy, uz, vx, vy, vz
  real(long) :: gamd(3,3), gamu(3,3), ovdfc(3)
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
!
        psifc  =  psif(irf,itf,ipf)
        alpfc  = alphf(irf,itf,ipf)
        bvxdfc = bvxdf(irf,itf,ipf)
        bvydfc = bvydf(irf,itf,ipf)
        bvzdfc = bvzdf(irf,itf,ipf)
        gamd(1,1) = hxxdf(irf,itf,ipf) + 1.0d0
        gamd(1,2) = hxydf(irf,itf,ipf)
        gamd(1,3) = hxzdf(irf,itf,ipf)
        gamd(2,2) = hyydf(irf,itf,ipf) + 1.0d0
        gamd(2,3) = hyzdf(irf,itf,ipf)
        gamd(3,3) = hzzdf(irf,itf,ipf) + 1.0d0
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)
        gamu(1,1) = hxxuf(irf,itf,ipf) + 1.0d0
        gamu(1,2) = hxyuf(irf,itf,ipf)
        gamu(1,3) = hxzuf(irf,itf,ipf)
        gamu(2,2) = hyyuf(irf,itf,ipf) + 1.0d0
        gamu(2,3) = hyzuf(irf,itf,ipf)
        gamu(3,3) = hzzuf(irf,itf,ipf) + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
        psifc4 = psifc**4
        psifc4inv = 1.0d0/psifc**4
        psifc8 = psifc**8
!
        emdfc = emd(irf,itf,ipf)
        if (irf.eq.nrf) then
          emdfc = 0.0d0
        end if
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, enefc)
        utfc  = utf(irf,itf,ipf)
        uxfc  = uxf(irf,itf,ipf)
        uyfc  = uyf(irf,itf,ipf)
        uzfc  = uzf(irf,itf,ipf)
!
        ut = utfc
        ux = uxfc
        uy = uyfc
        uz = uzfc
        vx = ux/ut
        vy = uy/ut
        vz = uz/ut
        ovdfc(1) = bvxdfc + gamd(1,1)*vx + gamd(1,2)*vy + gamd(1,3)*vz
        ovdfc(2) = bvydfc + gamd(2,1)*vx + gamd(2,2)*vy + gamd(2,3)*vz
        ovdfc(3) = bvzdfc + gamd(3,1)*vx + gamd(3,2)*vy + gamd(3,3)*vz
!
        rhoH(irf,itf,ipf) = hhfc*rhofc*(alpfc*utfc)**2 - prefc
!
        do ia = 1, 3
          jmd(irf,itf,ipf,ia) = hhfc*rhofc*alpfc*utfc**2*psifc4*ovdfc(ia)
          do ib = 1, 3
            smijd(irf,itf,ipf,ia,ib) &
            &                 = hhfc*rhofc*utfc**2*psifc8*ovdfc(ia)*ovdfc(ib) &
            &                 + prefc*psifc4*gamd(ia,ib)
          end do
        end do
!
        trsm(irf,itf,ipf) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            trsm(irf,itf,ipf) = trsm(irf,itf,ipf) &
            &                 + psifc4inv*gamu(ia,ib)*smijd(irf,itf,ipf,ia,ib)
          end do
        end do
!
      end do
    end do
  end do
!
end subroutine SEM_tensor
