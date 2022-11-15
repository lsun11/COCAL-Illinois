subroutine calc_admmom_asympto(cobj)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgout
  use coordinate_grav_r, only : rg
  use def_matter_parameter, only : radi
  use make_array_2d
  use make_array_3d
  use def_quantities, only : admmom_asymp
  use interface_source_admmom_asympto
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac8pi
  real(long) :: surf, rg_asympt
  real(long), pointer :: sousf(:,:), sousfv(:,:,:)
  integer    :: irg, ir, ii
  character(len=2), intent(in) :: cobj
!
  call alloc_array2d(sousf , 0, ntg, 0, npg)
  call alloc_array3d(sousfv, 0, ntg, 0, npg, 1, 3)
!      
  rg_asympt = 10.0**(dlog10(rgout)*2.0/3.0)
  do ir = 0, nrg
    irg = ir
    if (rg(ir).ge.rg_asympt) exit
  end do
!
  call source_admmom_asympto(sousfv,irg)
  do ii = 1, 3
    sousf(0:ntg,0:npg) = sousfv(0:ntg,0:npg,ii)
    call surf_int_grav_rg(sousf,surf,irg)
    fac8pi = 1.0d0/(8.0d0*pi)
    admmom_asymp(ii) = fac8pi*radi**2*surf
    if (cobj.eq.'bh') admmom_asymp(ii) = fac8pi*surf
  end do
!
!  write (6,'(a20,1p,e14.6)') ' ADM J      =       ', angmom_asymp
!
  deallocate(sousf)
  deallocate(sousfv)
end subroutine calc_admmom_asympto
