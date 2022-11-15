subroutine calc_gradvep_from_corot_id(potf,potxf,potyf,potzf)
  use phys_constant, only  :   long
  use grid_parameter
  use def_metric_on_SFC_CF
  use def_matter 
  use coordinate_grav_r, only : rg
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use def_matter_parameter, only : ome, ber
  use def_vector_phi, only : vec_phif
!  use def_matter, only : vep
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  implicit none
  real(long), pointer :: potf(:,:,:), potxf(:,:,:), potyf(:,:,:), potzf(:,:,:)
  real(long) ::  qq, hh, pre, rho0, ene, ovdfc(3), vphif(3), ut, psi4
  integer    :: irf, itf, ipf
!
  call interpo_gr2fl_metric_CF

  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        qq = emd(irf,itf,ipf)
        call peos_q2hprho(qq, hh, pre, rho0, ene)
        
        vphif(1) = vec_phif(irf,itf,ipf,1)
        vphif(2) = vec_phif(irf,itf,ipf,2)
        vphif(3) = vec_phif(irf,itf,ipf,3)
        ovdfc(1) = bvxdf(irf,itf,ipf) + ome*vphif(1)
        ovdfc(2) = bvydf(irf,itf,ipf) + ome*vphif(2)
        ovdfc(3) = bvzdf(irf,itf,ipf) + ome*vphif(3)
        psi4 = psif(irf,itf,ipf)**4
        ut   = hh/ber

        potxf(irf,itf,ipf) = ovdfc(1)*psi4*hh*ut
        potyf(irf,itf,ipf) = ovdfc(2)*psi4*hh*ut
        potzf(irf,itf,ipf) = ovdfc(3)*psi4*hh*ut
      end do
    end do
  end do
!
  write(6,*) " ", ome, ber
  write(6,*) "gradvep x component:", potxf(nrf,ntf/2,0), potxf(0,0,0), potxf(nrf,ntf/2,npf/2)
  write(6,*) "gradvep y component:", potyf(nrf,ntf/2,0), potyf(0,0,0), potyf(nrf,ntf/2,npf/2)

end subroutine calc_gradvep_from_corot_id
