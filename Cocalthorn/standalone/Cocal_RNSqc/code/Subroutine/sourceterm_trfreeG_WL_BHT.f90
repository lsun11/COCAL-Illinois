subroutine sourceterm_trfreeG_WL_BHT(souten)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrgin, nrg, ntg, npg
  use def_metric, only : psi, alph,  bvxd, bvyd, bvzd
  use def_metric_hij
  use def_matter, only : omeg, emdg, jomeg_int
  use def_matter_parameter, only : ber, radi, emdc
  use def_vector_phi, only : hvec_phig
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: souten(:,:,:,:)
  real(long) :: gamu(1:3,1:3), gamd(1:3,1:3)
  real(long) :: sab(1:3,1:3), ovd(1:3), vphig(1:3)
  real(long) :: zfac, emdgc, hhgc, pregc, rhogc, &
  &             san, utgc, omegc, jomeg_intgc, &
  &             psigc4, sgam, tfsab, trsab, ene, &
  &             bvxdgc, bvydgc, bvzdgc, psigc, alpgc, &
  &             hxxdgc,hxydgc,hxzdgc,hyydgc,hyzdgc,hzzdgc,&
  &             hxxugc,hxyugc,hxzugc,hyyugc,hyzugc,hzzugc
  integer :: irg, itg, ipg, ia, ib, ic
!
! --- source for hij is evaluated on mid points.
!
  souten = 0.0d0
  zfac = emdc*1.0d-14
!
  san = 1.0d0/3.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = nrgin+1, nrg
        call interpo_linear_type0(emdgc,emdg,irg,itg,ipg)
        call interpo_linear_type0(omegc,omeg,irg,itg,ipg)
        call interpo_linear_type0(jomeg_intgc,jomeg_int,irg,itg,ipg)

        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxdgc,bvxd,irg,itg,ipg)
        call interpo_linear_type0(bvydgc,bvyd,irg,itg,ipg)
        call interpo_linear_type0(bvzdgc,bvzd,irg,itg,ipg)
        psigc4 = psigc**4

        call interpo_linear_type0(hxxdgc,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hxydgc,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hxzdgc,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hyydgc,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hyzdgc,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hzzdgc,hzzd,irg,itg,ipg)

        call interpo_linear_type0(hxxugc,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hxyugc,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hxzugc,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hyyugc,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hyzugc,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hzzugc,hzzu,irg,itg,ipg)

        gamu(1,1) = hxxugc + 1.0d0
        gamu(1,2) = hxyugc
        gamu(1,3) = hxzugc
        gamu(2,2) = hyyugc + 1.0d0
        gamu(2,3) = hyzugc
        gamu(3,3) = hzzugc + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
!
        gamd(1,1) = hxxdgc + 1.0d0
        gamd(1,2) = hxydgc
        gamd(1,3) = hxzdgc
        gamd(2,2) = hyydgc + 1.0d0
        gamd(2,3) = hyzdgc
        gamd(3,3) = hzzdgc + 1.0d0
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)
!
        if (emdgc > zfac) then
          call peos_q2hprho(emdgc, hhgc, pregc, rhogc, ene)
          utgc  = hhgc/ber*exp(jomeg_intgc)
          vphig(1) = hvec_phig(irg,itg,ipg,1)
          vphig(2) = hvec_phig(irg,itg,ipg,2)
          vphig(3) = hvec_phig(irg,itg,ipg,3)
        
! --- CF and waveless contributions
          ovd(1) = bvxdgc + gamd(1,1)*omegc*vphig(1) &
          &               + gamd(1,2)*omegc*vphig(2) &
          &               + gamd(1,3)*omegc*vphig(3)
          ovd(2) = bvydgc + gamd(2,1)*omegc*vphig(1) &
          &               + gamd(2,2)*omegc*vphig(2) &
          &               + gamd(2,3)*omegc*vphig(3)
          ovd(3) = bvzdgc + gamd(3,1)*omegc*vphig(1) &
          &               + gamd(3,2)*omegc*vphig(2) &
          &               + gamd(3,3)*omegc*vphig(3)

          do ic = 1, 6
            ia = 1 + ic/4 + ic/6
            ib = ic - (ic/4)*2 - ic/6
            sab(ia,ib) = hhgc*rhogc*utgc**2*psigc4**2*ovd(ia)*ovd(ib)
          end do
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
            souten(irg,itg,ipg,ic) = 2.0d0*(- radi**2*8.0d0*pi*tfsab)
          end do
        end if
!
      end do
    end do
  end do
!
end subroutine sourceterm_trfreeG_WL_BHT
