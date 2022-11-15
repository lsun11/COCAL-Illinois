subroutine sourceterm_trfreeG_WL_EMF(souten)
  use phys_constant, only : pi
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : hrg
  use def_metric, only     : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_SEM_tensor_EMF, only : smijd_EMF, trsm_EMF
  use interface_interpo_linear_type0
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8) :: gamd(1:3,1:3), sab(1:3,1:3)
  real(8) :: psigc, hhxxd, hhxyd, hhxzd, hhyyd, hhyzd, hhzzd, &
  &          san = 1.0d0/3.0d0, sgam, tfsab, trsab
  integer :: ipg, irg, itg, ia, ib, ic
!
! --- source for hij is evaluated on grid points.
!
!
!!      open(15,file='test_vec_sou_x',status='unknown')
!!      open(16,file='test_vec_sou_y',status='unknown')
!!      open(17,file='test_vec_sou_z',status='unknown')
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
! --- EMF source
        do ib = 1, 3
          do ia = 1, 3
            sab(ia,ib) = smijd_EMF(irg,itg,ipg,ia,ib)
          end do
        end do
!
        trsab = trsm_EMF(irg,itg,ipg)
!
        do ic = 1, 6
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
          sgam = san*gamd(ia,ib)
          tfsab = sab(ia,ib) - psigc**4*sgam*trsab
          souten(irg,itg,ipg,ic) = 2.0d0*(-8.0d0*pi*tfsab)
!
!
!!      if (itg.eq.73.and.ipg.eq.16.and.ic.eq.1) then 
!!          write(15,'(1p,20e20.12)')  hrg(irg),  &
!!             &      souten(irg,itg,ipg,ic) , &
!!             &  tfsab,  sab(ia,ib), - psigc**4*sgam*trsab
!!      end if
!!      if (itg.eq.73.and.ipg.eq.16.and.ic.eq.4) then 
!!          write(16,'(1p,20e20.12)')  hrg(irg),  &
!!             &      souten(irg,itg,ipg,ic) , &
!!             &  tfsab,  sab(ia,ib), - psigc**4*sgam*trsab
!!      end if
!!      if (itg.eq.73.and.ipg.eq.16.and.ic.eq.6) then 
!!          write(17,'(1p,20e20.12)')  hrg(irg),  &
!!             &      souten(irg,itg,ipg,ic) , &
!!             &  tfsab,  sab(ia,ib), - psigc**4*sgam*trsab
!!      end if
!
        end do
!
      end do
    end do
  end do
!
!!  close(15)
!!  close(16)
!!  close(17)
end subroutine sourceterm_trfreeG_WL_EMF
