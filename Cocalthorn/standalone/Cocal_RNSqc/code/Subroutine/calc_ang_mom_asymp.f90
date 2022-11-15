subroutine calc_ang_mom_asymp(cobj)
!
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgout
  use coordinate_grav_r, only : hrg
  use def_matter_parameter, only : radi
  use def_quantities, only : angmom_asymp
  use make_array_2d
  use interface_source_ang_mom_asymp
  use interface_surf_int_grav
  implicit none
  real(long) :: fac8pi
  real(long)          :: surf, rg_asympt
  real(long), pointer :: sousf(:,:)
  integer             :: irg, ir
  character(len=2), intent(in) :: cobj
!
  call alloc_array2d(sousf,1,ntg,1,npg)
!
!!  irg = nrg - 4
  rg_asympt = 10.0**(dlog10(rgout)*2.0/3.0)
  do ir = 1, nrg
    irg = ir
    if (hrg(ir).ge.rg_asympt) exit
  end do
!
  call source_ang_mom_asymp(sousf,irg)
  call surf_int_grav(sousf,surf,irg)
  fac8pi = 1.0d0/(8.0d0*pi)
  angmom_asymp = fac8pi*radi**2*surf
  if (cobj.eq.'bh') angmom_asymp = fac8pi*surf
!
!  write (6,'(a20,1p,e14.6)') ' Angular momentum = ', angmom_asymp
!
  deallocate(sousf)
!
end subroutine calc_ang_mom_asymp
