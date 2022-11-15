subroutine calc_ang_mom_BBH_CF_inf
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg
  use make_array_1d
  use make_array_2d
  use def_quantities, only : angmom, angmom_asymp
  use interface_source_ang_mom_inf
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac8pi
  real(long) :: int_ang_mom_inf
  real(long),pointer :: sou_ang_mom_inf(:,:)
  integer :: mass_ir, irg
!
  call alloc_array2d(sou_ang_mom_inf, 1, ntg, 1, npg)
!
  fac8pi = 0.125d0/pi
!
  call calc_mass_ir(mass_ir)
!  write(6,*) 'mass_ir = ', mass_ir
!
  call source_ang_mom_inf(sou_ang_mom_inf,mass_ir)

  call surf_int_grav_rg(sou_ang_mom_inf, int_ang_mom_inf, mass_ir)

  angmom_asymp = fac8pi*int_ang_mom_inf
  angmom = angmom_asymp
!
!  write(6,*) '---------------------------------------------------------------------------'
!  write(6,'(a28,1p,e14.6)') 'Angular momentum (infinity)=', angmom_asymp

  deallocate(sou_ang_mom_inf)
end subroutine calc_ang_mom_BBH_CF_inf
