subroutine source_ang_mom_asymp(sousf,irg)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg
  use def_metric, only  :   tfkij, psi, alph, bvxd, bvyd, bvzd
  use coordinate_grav_r, only       : hrg
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only   : hsinphig, hcosphig
  use def_vector_phi, only : hvec_phig
  implicit none
  real(long), pointer :: sousf(:,:)
  integer, intent(in) :: irg
  integer             :: itg, ipg, ia, ib
  real(long)  ::   vphig(1:3), rna(1:3)
  real(long)  ::   psiw
!
  do ipg = 1, npg
    do itg = 1, ntg
!
      vphig(1) = hvec_phig(irg,itg,ipg,1)
      vphig(2) = hvec_phig(irg,itg,ipg,2)
      vphig(3) = hvec_phig(irg,itg,ipg,3)
      rna(1) = hsinthg(itg)*hcosphig(ipg)
      rna(2) = hsinthg(itg)*hsinphig(ipg)
      rna(3) = hcosphig(ipg)
!
      sousf(itg,ipg) = 0.0d0
      do ib = 1, 3
        do ia = 1, 3
          sousf(itg,ipg) = sousf(itg,ipg)  &
          &              + tfkij(irg,itg,ipg,ia,ib)*vphig(ia)*rna(ib)
        end do
      end do
!
    end do
  end do
!
end subroutine source_ang_mom_asymp
