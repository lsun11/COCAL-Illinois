subroutine SEM_tensor_EMF
  use phys_constant,  only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,     only : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_faraday_tensor, only : fxd, fyd, fzd, fijdu, fijd, fidfiu, fijfij
  use def_SEM_tensor_EMF, only : rhoH_EMF, jmd_EMF, smijd_EMF, trsm_EMF
  use interface_interpo_linear_type0
  use make_array_3d
  implicit none
  integer :: irg, itg, ipg, ia, ib
  real(long) :: pi4inv, pi8inv
  real(long) :: psigc, psigc4, psigc4inv, psigc8inv, gamd(3,3), gamu(3,3)
  real(long) :: hxxdgc, hxydgc, hxzdgc, hyydgc, hyzdgc, hzzdgc
  real(long) :: hxxugc, hxyugc, hxzugc, hyyugc, hyzugc, hzzugc
  real(long) :: fidfiugc, fijfijgc, fid(3), fijdugc(3,3), fijdgc(3,3)
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  pi4inv = 1.0d0/(4.0d0*pi)
  pi8inv = 1.0d0/(8.0d0*pi)
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psigc,psi ,irg,itg,ipg)
        call interpo_linear_type0(hxxdgc,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hxydgc,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hxzdgc,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hyydgc,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hyzdgc,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hzzdgc,hzzd,irg,itg,ipg)
        psigc4 = psigc**4
        psigc4inv = 1.0d0/psigc**4
        psigc8inv = 1.0d0/psigc**8
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
        fidfiugc = fidfiu(irg,itg,ipg)
        fijfijgc = fijfij(irg,itg,ipg)
        fid(1) = fxd(irg,itg,ipg)
        fid(2) = fyd(irg,itg,ipg)
        fid(3) = fzd(irg,itg,ipg)
        fijdugc(1:3,1:3) = fijdu(irg,itg,ipg,1:3,1:3)
        fijdgc(1,1) = 0.0d0 ; fijdgc(2,2) = 0.0d0 ; fijdgc(3,3) = 0.0d0
        fijdgc(1,2) = fijd(irg,itg,ipg,1) ; fijdgc(2,1) = - fijdgc(1,2)
        fijdgc(1,3) = fijd(irg,itg,ipg,2) ; fijdgc(3,1) = - fijdgc(1,3)
        fijdgc(2,3) = fijd(irg,itg,ipg,3) ; fijdgc(3,2) = - fijdgc(2,3)
!
        rhoH_EMF(irg,itg,ipg) = pi8inv*(psigc4inv*fidfiugc &
        &                     +  0.5d0* psigc8inv*fijfijgc)
!
        do ia = 1, 3
!
          jmd_EMF(irg,itg,ipg,ia) = pi4inv*psigc4inv &
          &                       *(fijdugc(ia,1)*fid(1) &
          &                       + fijdugc(ia,2)*fid(2) &
          &                       + fijdugc(ia,3)*fid(3))
!
          do ib = 1, 3
            smijd_EMF(irg,itg,ipg,ia,ib) = pi4inv*              &
            &         (  psigc4inv*(fijdugc(ia,1)*fijdgc(ib,1)  &
            &                     + fijdugc(ia,2)*fijdgc(ib,2)  &
            &                     + fijdugc(ia,3)*fijdgc(ib,3)) &
            &                     - fid(ia)*fid(ib)             &
            &          - 0.25d0*psigc4*gamd(ia,ib)              &
            &         *(psigc8inv*fijfijgc - 2.0d0*psigc4inv*fidfiugc) )
          end do
        end do
!
        trsm_EMF(irg,itg,ipg) = pi8inv*(psigc4inv*fidfiugc &
        &                     +  0.5d0* psigc8inv*fijfijgc)
!
      end do
    end do
  end do
!
! rhoH_EMF=0.0 ; jmd_EMF=0.0 ; smijd_EMF=0.0 ; trsm_EMF=0.0; 
!
end subroutine SEM_tensor_EMF
