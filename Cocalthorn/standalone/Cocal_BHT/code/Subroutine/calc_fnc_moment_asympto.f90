subroutine calc_fnc_moment_asympto(fnc,fnc_moment)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgout
  use coordinate_grav_r, only : rg
  use make_array_2d
  use make_array_3d
  use interface_source_fnc_moment_asympto
  use interface_surf_int_grav_solidangle
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long) :: fnc_moment(3)
  real(long) :: fac4pi
  real(long) :: surf, rg_asympt
  real(long), pointer :: sousf(:,:), sousfv(:,:,:)
  integer    :: irg, ir, ii
!
  call alloc_array2d(sousf , 0, ntg, 0, npg)
  call alloc_array3d(sousfv, 0, ntg, 0, npg, 1, 3)
!      
  rg_asympt = 10.0**(dlog10(rgout)*2.0/3.0)
!  rg_asympt = 100.0d0
  do ir = 0, nrg
    irg = ir
    if (rg(ir).ge.rg_asympt) exit
  end do
!
  call source_fnc_moment_asympto(fnc,sousfv,irg)
  do ii = 1, 3
    sousf(0:ntg,0:npg) = sousfv(0:ntg,0:npg,ii)
    call surf_int_grav_solidangle(sousf,surf)
    fac4pi = 1.0d0/(4.0d0*pi)
    fnc_moment(ii) = fac4pi*surf
  end do
!
!  write (6,'(a20,1p,e14.6)') ' ADM J      =       ', angmom_asymp
!
  deallocate(sousf)
  deallocate(sousfv)
end subroutine calc_fnc_moment_asympto
