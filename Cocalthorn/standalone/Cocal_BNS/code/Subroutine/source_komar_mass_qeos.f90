subroutine source_komar_mass_qeos(soug,souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   rhof, utf
  use def_matter_parameter, only  :   radi, ber, rhos_qs
  use def_metric, only  :   tfkijkij,  psi, alph
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_interpo_linear_type0
  implicit none
  real(long),pointer ::   soug(:,:,:), souf(:,:,:)
  real(long), pointer ::  alphf(:,:,:), psif(:,:,:) 
  real(long)  ::   psiwm6
  real(long)  ::    alphw, psiw, rhow, prew, hhw, utw, rhoHw, esseS
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   otermx, otermy, otermz, ene, dummy
  integer     ::   irg,itg,ipg,irf,itf,ipf
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
        call interpo_linear_type0(alphw,alph,irg,itg,ipg)
        psiwm6 = psiw**6
        soug(irg,itg,ipg) = alphw*psiwm6*tfkijkij(irg,itg,ipg)
      end do
    end do
  end do
!
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        rhow = rhof(irf,itf,ipf)
        psiw = psif(irf,itf,ipf)
        alphw = alphf(irf,itf,ipf)
        call quark_rho2phenedpdrho(rhow, prew, hhw, ene, dummy)
        utw = utf(irf,itf,ipf)
!
        rhoHw = hhw*rhow*(alphw*utw)**2 - prew
        esseS = -hhw*rhow + 4.0d0*prew + rhoHw
! 
        souf(irf,itf,ipf) = 4.0d0*pi*alphw*psiw**6*(esseS+rhoHw)
      end do
    end do
  end do
!
  deallocate(alphf)
  deallocate(psif)
end subroutine source_komar_mass_qeos
