subroutine source_proper_mass(souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   emdg, emd, utf
  use def_matter_parameter, only  :   radi, pinx, ber
  use def_metric, only  :   tfkijkij,  psi, alph
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long),pointer ::   souf(:,:,:)
  real(long), pointer ::  alphf(:,:,:), psif(:,:,:) 
  real(long)  ::   psiwm6
  real(long)  ::   emdw, alphw, psiw, epsilonw, hhw, utw, rhoHw, esseS
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   otermx, otermy, otermz
  integer     ::   ir,it,ip
!
  call alloc_array3d(psif, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(alphf, 0, nrf, 0, ntf, 0, npf)
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
!
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        emdw = emd(ir,it,ip)
        if (emdw <= small) emdw = small
        psiw = psif(ir,it,ip)
        alphw = alphf(ir,it,ip)
        epsilonw = emdw**pinx*(1.0d0+pinx*emdw)
        hhw  = 1.0d0 + (pinx+1.0d0)*emdw
!        utw = hhw/ber
        utw = utf(irf,itf,ipf)
!
        souf(ir,it,ip) = epsilonw*alphw*utw*psiw**6
      end do
    end do
  end do
!
  deallocate(alphf)
  deallocate(psif)
end subroutine source_proper_mass
