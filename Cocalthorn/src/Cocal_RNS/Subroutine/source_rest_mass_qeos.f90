subroutine source_rest_mass_qeos(souf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only  :   rhof, utf
  use def_matter_parameter, only  :   radi, ber, rhos_qs
  use def_metric, only  :   tfkijkij,  psi, alph
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long),pointer ::   souf(:,:,:)
  real(long), pointer ::  alphf(:,:,:), psif(:,:,:) 
  real(long)  ::   psiwm6
  real(long)  ::    alphw, psiw, rhow,  utw, rhoHw, esseS
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
        rhow = rhof(ir,it,ip)
        if (rhow <= rhos_qs) rhow = small
        psiw = psif(ir,it,ip)
        alphw = alphf(ir,it,ip)
        utw = utf(ir,it,ip)
!
        souf(ir,it,ip) = rhow*alphw*utw*psiw**6
      end do
    end do
  end do
!
  deallocate(alphf)
  deallocate(psif)
end subroutine source_rest_mass_qeos
