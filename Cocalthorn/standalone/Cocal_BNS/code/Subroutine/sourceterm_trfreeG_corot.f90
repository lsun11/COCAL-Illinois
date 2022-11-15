subroutine sourceterm_trfreeG_corot(souten)
  use phys_constant, only : pi
  use def_matter_parameter, only : ber, radi
  use def_metric, only : psi, alph
  use def_metric_rotshift, only : ovxd, ovyd, ovzd
  use def_matter, only : emdg
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, &
  &                          hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use interface_interpo_linear_type0
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8) :: gamu(1:3,1:3), gamd(1:3,1:3)
  real(8) :: sab(1:3,1:3), ovd(1:3)
!
  real(8) :: emdgc, hhgc, pregc, psigc, rhogc, &
  &          hhxxu, hhxyu, hhxzu, hhyyu, hhyzu, hhzzu, &
  &          san, utgc, zfac, &
  &          hhxxd, hhxyd, hhxzd, hhyyd, hhyzd, hhzzd, &
  &          psigc4, sgam, tfsab, trsab, ene, &
  &          ovxgc, ovygc, ovzgc
  integer :: ipg, irg, itg, ia, ib, ic
!
! --- source for hij is evaluated on grid points.
!
  san = 1.0d0/3.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        psigc4 = psigc**4
        call interpo_linear_type0(emdgc,emdg,irg,itg,ipg)
        call peos_q2hprho(emdgc,hhgc,pregc,rhogc,ene)
        utgc  = hhgc/ber
!
        zfac = 1.0d0
        if (emdgc <= 1.0d-14) zfac = 0.0d0
!
        call interpo_linear_type0(hhxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hhxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hhxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hhyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hhyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hhzzu,hzzu,irg,itg,ipg)
        gamu(1,1) = hhxxu + 1.0d0
        gamu(1,2) = hhxyu
        gamu(1,3) = hhxzu
        gamu(2,2) = hhyyu + 1.0d0
        gamu(2,3) = hhyzu
        gamu(3,3) = hhzzu + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
!
        call interpo_linear_type0(hhxxd,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hhxyd,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hhxzd,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hhyyd,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hhyzd,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hhzzd,hzzd,irg,itg,ipg)
        gamd(1,1) = hhxxd + 1.0d0
        gamd(1,2) = hhxyd
        gamd(1,3) = hhxzd
        gamd(2,2) = hhyyd + 1.0d0
        gamd(2,3) = hhyzd
        gamd(3,3) = hhzzd + 1.0d0
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)
!
        call interpo_linear_type0(ovxgc,ovxd,irg,itg,ipg)
        call interpo_linear_type0(ovygc,ovyd,irg,itg,ipg)
        call interpo_linear_type0(ovzgc,ovzd,irg,itg,ipg)
        ovd(1) = ovxgc
        ovd(2) = ovygc
        ovd(3) = ovzgc
!
        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
! --- Fluid source
          sab(ia,ib) = hhgc*rhogc*utgc**2*psigc4**2*ovd(ia)*ovd(ib)
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
          souten(irg,itg,ipg,ic) = 2.0d0*(- radi**2*8.0d0*pi*tfsab*zfac)
        end do
      end do
    end do
  end do
!
end subroutine sourceterm_trfreeG_corot
