subroutine source_adm_mass_peos(soug,souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, utf
  use def_matter_parameter, only  :  radi, ber
  use def_metric, only  :   tfkijkij, psi, alph
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: soug(:,:,:)
  real(long), pointer :: souf(:,:,:)
  real(long), pointer :: alphf(:,:,:), psif(:,:,:) 
  integer     ::   irg,itg,ipg,irf,itf,ipf
  real(long)  ::   psiwm7
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, rhoHw
  real(long)  ::   epsilonw, alutw
  real(long)  ::   zfac, small = 1.0d-15
!
  call alloc_array3d(psif, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf)
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psiw,psi,irg,itg,ipg)
        soug(irg,itg,ipg) = 0.125d0*psiw**5*tfkijkij(irg,itg,ipg)
      end do
    end do
  end do
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        emdw = emd(irf,itf,ipf)
        if (emdw <= small) emdw = small
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        call peos_q2hprho(emdw, hhw, prew, rhow, epsilonw)
        utw = utf(irf,itf,ipf)
!        rhoHw = hhw*rhow*(alphw*utw)**2 - prew
        alutw = alphw*utw
!
!        souf(irf,itf,ipf) = 2.0d0*pi*psiw**5*rhoHw
        souf(irf,itf,ipf) = 2.0d0*pi*psiw**5  &
        &  *(epsilonw*alutw**2 + prew*(alutw-1.0d0)*(alutw+1.0d0))
      end do
    end do
  end do
!
  deallocate(alphf)
  deallocate(psif)
end subroutine source_adm_mass_peos
