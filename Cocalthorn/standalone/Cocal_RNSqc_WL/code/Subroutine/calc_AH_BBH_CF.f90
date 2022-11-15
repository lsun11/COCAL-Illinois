subroutine calc_AH_BBH_CF
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   ntg, npg
  use def_quantities_bh, only : AHmass, AHarea
  use interface_source_AHarea_CF
  use interface_surf_int_grav_rg
  use make_array_2d
  implicit none
  real(long) :: surf
  real(long), pointer :: sousf(:,:)
!
  call alloc_array2d(sousf,0,ntg,0,npg)
!
  call source_AHarea_CF(sousf)
  call surf_int_grav_rg(sousf,surf,0)
  AHarea = surf
  AHmass = dsqrt(AHarea/(16.0d0*pi))
!
  deallocate(sousf)
end subroutine calc_AH_BBH_CF
