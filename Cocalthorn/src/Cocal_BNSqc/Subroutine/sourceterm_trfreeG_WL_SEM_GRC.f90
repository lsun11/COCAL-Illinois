subroutine sourceterm_trfreeG_WL_SEM(souten)
  use phys_constant, only : pi
  use def_matter_parameter, only : radi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric    , only : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_SEM_tensor, only : smijd, trsm
  use interface_interpo_linear_type0
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8) :: gamd(1:3,1:3), sab(1:3,1:3)
  real(8) :: psigc, hhxxd, hhxyd, hhxzd, hhyyd, hhyzd, hhzzd, &
  &          san, sgam, tfsab, trsab
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
        call interpo_linear_type0(psigc,psi ,irg,itg,ipg)
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
! --- Fluid source
        do ib = 1, 3
          do ia = 1, 3
            sab(ia,ib) = smijd(irg,itg,ipg,ia,ib)
          end do
        end do
!
        trsab = trsm(irg,itg,ipg)
!
        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
          sgam = san*gamd(ia,ib)
          tfsab = sab(ia,ib) - psigc**4*sgam*trsab
          souten(irg,itg,ipg,ic) = 2.0d0*(- radi**2*8.0d0*pi*tfsab)
        end do
      end do
    end do
  end do
!
end subroutine sourceterm_trfreeG_WL_SEM
