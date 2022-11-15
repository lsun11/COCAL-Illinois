subroutine sourceterm_trfreeG_drot_SFC(souten)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF, only : psif, bvxdf, bvydf, bvzdf
  use def_metric_on_SFC_WL, only : hxxdf, hxydf, hxzdf, hyydf, hyzdf, hzzdf, &
  &                                hxxuf, hxyuf, hxzuf, hyyuf, hyzuf, hzzuf
!  use def_metric_rotshift, only : ovxd, ovyd, ovzd
  use def_matter, only : emd, omef, jomef_int
  use def_matter_parameter, only : ber, radi
  use def_vector_phi, only : vec_phif
  implicit none
  real(long), pointer :: souten(:,:,:,:)
  real(long) :: gamu(1:3,1:3), gamd(1:3,1:3)
  real(long) :: sab(1:3,1:3), ovd(1:3), vphif(1:3)
!
  real(long) :: emdfc, hhfc, prefc, psifc, rhofc, &
  &             san, utfc, omefc, jomef_intfc, &
  &             psifc4, sgam, tfsab, trsab, ene, &
  &             bvxfc, bvyfc, bvzfc
  integer :: irf, itf, ipf, ia, ib, ic
!
! --- source for hij is evaluated on grid points.
!
  san = 1.0d0/3.0d0
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
!
        gamu(1,1) = hxxuf(irf,itf,ipf) + 1.0d0
        gamu(1,2) = hxyuf(irf,itf,ipf)
        gamu(1,3) = hxzuf(irf,itf,ipf)
        gamu(2,2) = hyyuf(irf,itf,ipf) + 1.0d0
        gamu(2,3) = hyzuf(irf,itf,ipf)
        gamu(3,3) = hzzuf(irf,itf,ipf) + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
!
        gamd(1,1) = hxxdf(irf,itf,ipf) + 1.0d0
        gamd(1,2) = hxydf(irf,itf,ipf)
        gamd(1,3) = hxzdf(irf,itf,ipf)
        gamd(2,2) = hyydf(irf,itf,ipf) + 1.0d0
        gamd(2,3) = hyzdf(irf,itf,ipf)
        gamd(3,3) = hzzdf(irf,itf,ipf) + 1.0d0
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)
!
        bvxfc = bvxdf(irf,itf,ipf)
        bvyfc = bvydf(irf,itf,ipf)
        bvzfc = bvzdf(irf,itf,ipf)
!
        psifc  = psif(irf,itf,ipf)
        psifc4 = psifc**4
        emdfc  =  emd(irf,itf,ipf)
        omefc  = omef(irf,itf,ipf)
        jomef_intfc = jomef_int(irf,itf,ipf)
        if (irf.eq.nrf) then
          emdfc = 0.0d0
        end if
        call peos_q2hprho(emdfc,hhfc,prefc,rhofc,ene)
        utfc  = hhfc/ber*exp(jomef_intfc)
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
!
        ovd(1) = bvxfc + gamd(1,1)*omefc*vphif(1) &
        &              + gamd(1,2)*omefc*vphif(2) &
        &              + gamd(1,3)*omefc*vphif(3)
        ovd(2) = bvyfc + gamd(2,1)*omefc*vphif(1) &
        &              + gamd(2,2)*omefc*vphif(2) &
        &              + gamd(2,3)*omefc*vphif(3)
        ovd(3) = bvzfc + gamd(3,1)*omefc*vphif(1) &
        &              + gamd(3,2)*omefc*vphif(2) &
        &              + gamd(3,3)*omefc*vphif(3)
!
        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
! --- Fluid source
          sab(ia,ib) = hhfc*rhofc*utfc**2*psifc4**2*ovd(ia)*ovd(ib)
        end do
! --- end
        sab(2,1) = sab(1,2)
        sab(3,1) = sab(1,3)
        sab(3,2) = sab(2,3)
!
        trsab = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            trsab = trsab + gamu(ia,ib) * sab(ia,ib)
          end do
        end do
!
        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
          sgam = san*gamd(ia,ib)
          tfsab = sab(ia,ib) - sgam*trsab
          souten(irf,itf,ipf,ic) = 2.0d0*(- radi**2*8.0d0*pi*tfsab)
        end do
      end do
    end do
  end do
!
end subroutine sourceterm_trfreeG_drot_SFC
