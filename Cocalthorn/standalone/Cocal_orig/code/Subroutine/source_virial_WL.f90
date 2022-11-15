subroutine source_virial_WL(sou_Tkin,sou_Pint,sou_Wgra)
  use phys_constant,  only : long
  use interface_source_virial_matter
  use interface_source_virial_gravity_WL
  implicit none
  real(long), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
!
! --- Compute source terms for virial relations.
! --- sources for kinetic and internal energies are defined on SCF grid points.
! --- sources for EM, and gravitational energies are defined on mid points.
!
  call source_virial_matter(sou_Tkin,sou_Pint)
  call source_virial_gravity_WL(sou_Wgra)
!
end subroutine source_virial_WL
