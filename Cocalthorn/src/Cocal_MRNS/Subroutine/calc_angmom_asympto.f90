subroutine calc_angmom_asympto(cobj)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgout
  use coordinate_grav_r, only : rg
  use def_matter_parameter, only : radi
  use make_array_2d
  use def_quantities, only : angmom_asymp
  use interface_source_angmom_asympto
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac8pi
  real(long) :: surf, rg_asympt
  real(long), pointer :: sousf(:,:)
  integer    :: irg, ir
  character(len=2), intent(in) :: cobj
!
  call alloc_array2d(sousf, 0, ntg, 0, npg)
!      
  rg_asympt = 10.0**(dlog10(rgout)*2.0/3.0)
  do ir = 0, nrg
    irg = ir
    if (rg(ir).ge.rg_asympt) exit
  end do
!
  call source_angmom_asympto(sousf,irg)
  call surf_int_grav_rg(sousf,surf,irg)
  fac8pi = 1.0d0/(8.0d0*pi)
  angmom_asymp = fac8pi*radi**2*surf
  if (cobj.eq.'bh') angmom_asymp = fac8pi*surf
!
!  write (6,'(a20,1p,e14.6)') ' ADM J      =       ', angmom_asymp
!
  deallocate(sousf)
end subroutine calc_angmom_asympto
