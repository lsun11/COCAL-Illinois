subroutine SEM_tensor_GRC
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,     only : psi, alph, bvxd, bvyd, bvzd
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_matter,     only : emd , utf, uxf, uyf, uzf, &
  &                          emdg, utg, uxg, uyg, uzg
  use def_SEM_tensor, only : rhoH, jmd, smijd, trsm
  use interface_grgrad_midpoint
  use interface_interpo_linear_type0
  use interface_interpo_fl2gr
  use make_array_3d
  implicit none
  integer :: irg, itg, ipg, ia, ib
  real(long) :: psigc, psigc4, psigc4inv, psigc8, &
  &             alpgc, bvxdgc, bvydgc, bvzdgc, &
  &             hxxdgc, hxydgc, hxzdgc, hyydgc, hyzdgc, hzzdgc, &
  &             hxxugc, hxyugc, hxzugc, hyyugc, hyzugc, hzzugc
  real(long) :: emdgc, hhgc, pregc, rhogc, enegc
  real(long) :: utgc, uxgc, uygc, uzgc, ut, ux, uy, uz, vx, vy, vz
  real(long) :: gamd(3,3), gamu(3,3), ovdgc(3)
  real(long), pointer :: zfac(:,:,:)
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  call alloc_array3d(zfac, 0, nrg, 0, ntg, 0, npg)
!
  call interpo_fl2gr(emd, emdg)
  call interpo_fl2gr(utf, utg)
  call interpo_fl2gr(uxf, uxg)
  call interpo_fl2gr(uyf, uyg)
  call interpo_fl2gr(uzf, uzg)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(emdgc,emdg,irg,itg,ipg)
        zfac(irg,itg,ipg) = 1.0d0
        if (emdgc <= 1.0d-15) then
          emdgc = 1.0d-15
          zfac(irg,itg,ipg) = 0.0d0
          cycle
        end if
        call peos_q2hprho(emdgc, hhgc, pregc, rhogc, enegc)
!
        call interpo_linear_type0(utgc ,utg ,irg,itg,ipg)
        call interpo_linear_type0(uxgc ,uxg ,irg,itg,ipg)
        call interpo_linear_type0(uygc ,uyg ,irg,itg,ipg)
        call interpo_linear_type0(uzgc ,uzg ,irg,itg,ipg)
        call interpo_linear_type0(psigc,psi ,irg,itg,ipg)
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxdgc,bvxd,irg,itg,ipg)
        call interpo_linear_type0(bvydgc,bvyd,irg,itg,ipg)
        call interpo_linear_type0(bvzdgc,bvzd,irg,itg,ipg)
        call interpo_linear_type0(hxxdgc,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hxydgc,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hxzdgc,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hyydgc,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hyzdgc,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hzzdgc,hzzd,irg,itg,ipg)
        psigc4 = psigc**4
        psigc4inv = 1.0d0/psigc**4
        psigc8 = psigc**8
        gamd(1,1) = 1.0d0 + hxxdgc
        gamd(1,2) =         hxydgc
        gamd(1,3) =         hxzdgc
        gamd(2,2) = 1.0d0 + hyydgc
        gamd(2,3) =         hyzdgc
        gamd(3,3) = 1.0d0 + hzzdgc
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)
!
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
        ut = utgc
        ux = uxgc
        uy = uygc
        uz = uzgc
        vx = ux/ut
        vy = uy/ut
        vz = uz/ut
        ovdgc(1) = bvxdgc + gamd(1,1)*vx + gamd(1,2)*vy + gamd(1,3)*vz
        ovdgc(2) = bvydgc + gamd(2,1)*vx + gamd(2,2)*vy + gamd(2,3)*vz
        ovdgc(3) = bvzdgc + gamd(3,1)*vx + gamd(3,2)*vy + gamd(3,3)*vz
!
        rhoH(irg,itg,ipg) = hhgc*rhogc*(alpgc*utgc)**2 - pregc
!
        do ia = 1, 3
          jmd(irg,itg,ipg,ia) = hhgc*rhogc*alpgc*utgc**2*psigc4*ovdgc(ia)
          do ib = 1, 3
            smijd(irg,itg,ipg,ia,ib) &
            &                 = hhgc*rhogc*utgc**2*psigc8*ovdgc(ia)*ovdgc(ib) &
            &                 + pregc*psigc4*gamd(ia,ib)
          end do
        end do
!
        trsm(irg,itg,ipg) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            trsm(irg,itg,ipg) = trsm(irg,itg,ipg) &
            &                 + psigc4inv*gamu(ia,ib)*smijd(irg,itg,ipg,ia,ib)
          end do
        end do
!
      end do
    end do
  end do
!
  rhoH(1:nrg,1:ntg,1:npg) = rhoH(1:nrg,1:ntg,1:npg) &
  &                        *zfac(1:nrg,1:ntg,1:npg)
  trsm(1:nrg,1:ntg,1:npg) = trsm(1:nrg,1:ntg,1:npg) &
  &                        *zfac(1:nrg,1:ntg,1:npg)
  do ia = 1, 3
    jmd(1:nrg,1:ntg,1:npg,ia) = jmd(1:nrg,1:ntg,1:npg,ia) &
    &                         *zfac(1:nrg,1:ntg,1:npg)
    do ib = 1, 3
      smijd(1:nrg,1:ntg,1:npg,ia,ib) = smijd(1:nrg,1:ntg,1:npg,ia,ib) &
      &                                *zfac(1:nrg,1:ntg,1:npg)
    end do
  end do
!
  deallocate(zfac)
!
end subroutine SEM_tensor_GRC
