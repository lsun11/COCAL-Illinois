subroutine calc_mass_asympto(cobj)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgout
  use coordinate_grav_r, only : rg
  use def_metric, only  :   psi, alph, alps
  use def_matter_parameter, only : radi
  use make_array_2d
  use def_quantities, only : admmass_asymp, komarmass_asymp
  use interface_source_mass_asympto
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
  call source_mass_asympto(psi,sousf,irg)
  call surf_int_grav_rg(sousf,surf,irg)
  fac2pi = 0.5d0/pi
  admmass_asymp = - fac2pi*radi*surf
  if (cobj.eq.'bh') admmass_asymp = - fac2pi*surf
!
  call source_mass_asympto(alph,sousf,irg)
  call surf_int_grav_rg(sousf,surf,irg)
  fac4pi = 0.25d0/pi
  komarmass_asymp = fac4pi*radi*surf
  if (cobj.eq.'bh') komarmass_asymp = fac4pi*surf
!
!  call source_mass_asympto(alps,sousf,irg)
!  call surf_int_grav_rg(sousf,surf,irg)
!
!  write (6,'(a20,1p,e14.6)') ' ADM   mass =       ', admmass
!  write (6,'(a20,1p,e14.6)') ' Komar mass =       ', komarmass
!
  deallocate(sousf)
end subroutine calc_mass_asympto
