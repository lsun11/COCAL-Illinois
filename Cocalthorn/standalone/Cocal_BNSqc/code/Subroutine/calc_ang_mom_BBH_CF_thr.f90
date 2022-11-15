subroutine calc_ang_mom_BBH_CF_thr
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg
  use make_array_1d
  use make_array_2d
  use def_quantities, only : angmom_thr
  use interface_source_ang_mom_thr
  use interface_source_ang_mom_exc
  use interface_surf_int_grav_rg
  use grid_parameter_binary_excision, only: ex_nrg
  implicit none
  real(long) :: fac8pi
  real(long) :: int_ang_mom_thr, int_ang_mom_exc
  real(long),pointer :: sou_ang_mom_thr(:,:), sou_ang_mom_exc(:,:)
  integer ::  irg
!
  call alloc_array2d(sou_ang_mom_thr, 1, ntg, 1, npg)
  call alloc_array2d(sou_ang_mom_exc, 1, ntg, 1, npg)
!
  fac8pi = 0.125d0/pi

  call source_ang_mom_thr(sou_ang_mom_thr)

  call surf_int_grav_rg(sou_ang_mom_thr, int_ang_mom_thr, 0)
!
! the following 3 lines are probably not needed...............................
  call source_ang_mom_exc(sou_ang_mom_exc)
  call surf_int_grav_rg(sou_ang_mom_exc, int_ang_mom_exc, ex_nrg)
  angmom_thr = -fac8pi*(int_ang_mom_thr + int_ang_mom_exc)
  write(6,'(a37,1p,e14.6)') 'Angular momentum (throat + excised) =', angmom_thr
!
  angmom_thr = -2.0d0*fac8pi*int_ang_mom_thr
  write(6,'(a28,1p,e14.6)') 'Angular momentum (throats) =', angmom_thr

  deallocate(sou_ang_mom_thr)
  deallocate(sou_ang_mom_exc)
end subroutine calc_ang_mom_BBH_CF_thr
