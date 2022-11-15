subroutine source_virial_matter(sou_Tkin,sou_Pint)
  use phys_constant,  only : long, pi
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF, only : psif, alphf
  use def_matter,     only : emd, utf
  use make_array_3d
  use make_array_4d
  implicit none
  integer :: irf, itf, ipf
  real(long) :: psifc, psifc6, alpfc, pi4inv
  real(long) :: emdfc, hhfc, prefc, rhofc, enefc, utfc
  real(long), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:)
!
! --- Compute matter source terms for virial relations.
! --- sources for kinetic and internal energies are defined on SCF grid points.
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
!
        emdfc = emd(irf,itf,ipf)
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, enefc)
        utfc = utf(irf,itf,ipf)
        psifc = psif(irf,itf,ipf)
        alpfc = alphf(irf,itf,ipf)
        psifc6 = psifc**6
!
        sou_Tkin(irf,itf,ipf)= 0.5d0*hhfc*rhofc*((alpfc*utfc)**2-1.0d0)*psifc6
        sou_Pint(irf,itf,ipf)= prefc*psifc6
!
      end do
    end do
  end do
!
end subroutine source_virial_matter
