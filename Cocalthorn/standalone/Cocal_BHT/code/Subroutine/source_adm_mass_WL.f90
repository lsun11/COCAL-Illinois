subroutine source_adm_mass_WL(soug,souf)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only     : emdg, emd, utf
  use def_matter_parameter, only : radi, ber
  use def_metric, only     : psi, alph, tfkijkij
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_metric_on_SFC_CF, only : psif, alphf
  use def_ricci_tensor, only : rab
  use def_SEM_tensor, only : rhoH
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: soug(:,:,:)
  real(long), pointer :: souf(:,:,:)
  integer     ::   irg, itg, ipg, irf, itf, ipf, ia, ib
  real(long)  ::   psiwm7
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, rhoHw
  real(long)  ::   epsilonw, alutw
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   rics, ric(1:3,1:3), gamu(1:3,1:3), &
  &                hxxuc, hxyuc, hxzuc, hyyuc, hyzuc, hzzuc
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        call interpo_linear_type0(hxxuc,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hxyuc,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hxzuc,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hyyuc,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hyzuc,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hzzuc,hzzu,irg,itg,ipg)
        gamu(1,1) = hxxuc + 1.0d0
        gamu(1,2) = hxyuc
        gamu(1,3) = hxzuc
        gamu(2,2) = hyyuc + 1.0d0
        gamu(2,3) = hyzuc
        gamu(3,3) = hzzuc + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
        ric(1,1) = rab(irg,itg,ipg,1)
        ric(1,2) = rab(irg,itg,ipg,2)
        ric(1,3) = rab(irg,itg,ipg,3)
        ric(2,2) = rab(irg,itg,ipg,4)
        ric(2,3) = rab(irg,itg,ipg,5)
        ric(3,3) = rab(irg,itg,ipg,6)
        ric(2,1) = ric(1,2)
        ric(3,1) = ric(3,1)
        ric(3,2) = ric(2,3)
        rics = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            rics = rics + gamu(ia,ib)*ric(ia,ib)
          end do
        end do
        soug(irg,itg,ipg) = - 0.125d0*psiw*rics &
        &                   + 0.125d0*psiw**5*tfkijkij(irg,itg,ipg)
      end do
    end do
  end do
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        psiw  = psif(irf,itf,ipf)
        rhoHw = rhoH(irf,itf,ipf)
        souf(irf,itf,ipf) = 2.0d0*pi*psiw**5*rhoHw
!!         emdw = emd(irf,itf,ipf)
!!         if (emdw <= small) emdw = small
!!         psiw = psif(irf,itf,ipf)
!!         alphw = alphf(irf,itf,ipf)
!!         call peos_q2hprho(emdw, hhw, prew, rhow, epsilonw)
!! !        utw = hhw/ber
!!         utw = utf(irf,itf,ipf)
!! !        rhoHw = hhw*rhow*(alphw*utw)**2 - prew
!! !        rhoHw = rhoH(irf,itf,ipf)
!!         alutw = alphw*utw
!! !
!! !        souf(irf,itf,ipf) = 2.0d0*pi*psiw**5*rhoHw
!!         souf(irf,itf,ipf) = 2.0d0*pi*psiw**5  &
!!         &  *(epsilonw*alutw**2 + prew*(alutw-1.0d0)*(alutw+1.0d0))
      end do
    end do
  end do
!
end subroutine source_adm_mass_WL
