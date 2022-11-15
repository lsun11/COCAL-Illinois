subroutine source_virial_CF_qeos(sou_Tkin,sou_Pint,sou_Wgra)
  use phys_constant,  only : long
  use interface_source_virial_matter_qeos
  use interface_source_virial_gravity_CF
  implicit none
  real(long), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
!
! --- Compute source terms for virial relations.
! --- sources for kinetic and internal energies are defined on SCF grid points.
! --- source  for gravitational energies is defined on mid points.
!
  call source_virial_matter_qeos(sou_Tkin,sou_Pint)
  call source_virial_gravity_CF(sou_Wgra)
!
end subroutine source_virial_CF_qeos
