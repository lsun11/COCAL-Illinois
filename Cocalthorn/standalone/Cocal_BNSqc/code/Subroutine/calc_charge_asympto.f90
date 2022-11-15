subroutine calc_charge_asympto(cobj)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgout
  use coordinate_grav_r, only : rg
  use def_matter_parameter, only : radi
  use make_array_2d
  use def_quantities, only : charge_asymp
  use interface_source_charge_asympto
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac2pi, fac4pi
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
  call source_charge_asympto(sousf,irg)
  call surf_int_grav_rg(sousf,surf,irg)
  fac4pi = 0.25d0/pi
  charge_asymp = fac4pi*radi*surf
  if (cobj.eq.'bh') charge_asymp = fac4pi*surf
!
  write (6,'(a20,1p,e14.6)') ' Charge asympt =    ', charge_asymp
!
  deallocate(sousf)
end subroutine calc_charge_asympto
