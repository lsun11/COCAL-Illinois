subroutine source_virial_CF(sou_Tkin,sou_Pint,sou_Wgra)
  use phys_constant,  only : long
  use interface_source_virial_matter
  use interface_source_virial_gravity_CF
  implicit none
  real(long), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
!
! --- Compute source terms for virial relations.
! --- sources for kinetic and internal energies are defined on SCF grid points.
! --- source  for gravitational energies is defined on mid points.
!
  call source_virial_matter(sou_Tkin,sou_Pint)
  call source_virial_gravity_CF(sou_Wgra)
!
end subroutine source_virial_CF
