subroutine sourceterm_trfreeG_WL_SEM(souten)
  use phys_constant, only : pi
  use def_matter_parameter, only : radi
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF, only : psif
  use def_metric_on_SFC_WL, only : hxxdf, hxydf, hxzdf, hyydf, hyzdf, hzzdf
  use def_SEM_tensor, only : smijd, trsm
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8) :: gamd(1:3,1:3), sab(1:3,1:3)
  real(8) :: psifc, san, sgam, tfsab, trsab
  integer :: ipf, irf, itf, ia, ib, ic
!
! --- source for hij is evaluated on grid points.
!
  san = 1.0d0/3.0d0
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
!
        psifc = psif(irf,itf,ipf)
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
! --- Fluid source
        do ib = 1, 3
          do ia = 1, 3
            sab(ia,ib) = smijd(irf,itf,ipf,ia,ib)
          end do
        end do
!
        trsab = trsm(irf,itf,ipf)
!
        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
          sgam = san*gamd(ia,ib)
          tfsab = sab(ia,ib) - psifc**4*sgam*trsab
          souten(irf,itf,ipf,ic) = 2.0d0*(- radi**2*8.0d0*pi*tfsab)
        end do
      end do
    end do
  end do
!
end subroutine sourceterm_trfreeG_WL_SEM
