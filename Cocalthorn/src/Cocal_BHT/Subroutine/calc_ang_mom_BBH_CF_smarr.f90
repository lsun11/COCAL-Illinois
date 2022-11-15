subroutine calc_ang_mom_BBH_CF_smarr
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg
  use make_array_1d
  use make_array_2d
  use def_quantities, only : admmass, angmom_smarr
  use def_bh_parameter, only : ome_bh
  use interface_source_ang_mom_smarr
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac4pi
  real(long) :: int_ang_mom_smarr, rhs
  real(long),pointer :: sou_ang_mom_smarr(:,:)
  integer ::  irg
!
  call alloc_array2d(sou_ang_mom_smarr, 1, ntg, 1, npg)
!
  fac4pi = 0.25d0/pi

  call source_ang_mom_smarr(sou_ang_mom_smarr)

  call surf_int_grav_rg(sou_ang_mom_smarr, int_ang_mom_smarr, 0)
 
  rhs = -2.0d0*fac4pi*int_ang_mom_smarr
!
!  write(6,'(a28,1p,e14.6)') 'RHS of Smarr formula       = ',rhs 

!  write(6,'(a10,1p,e14.6,a20,1p,e14.20)') 'admmass=', admmass, '            ome_bh=',ome_bh 
  angmom_smarr = (admmass - rhs)/2.0d0/ome_bh
!
!  write(6,'(a28,1p,e14.6,a10,a5,1p,e14.6,a1)') 'Angular momentum (Smarr)   =', angmom_smarr,  &
!        &      '', '(rhs=', rhs,')'
!  write(6,*) '---------------------------------------------------------------------------'

!

  deallocate(sou_ang_mom_smarr)
end subroutine calc_ang_mom_BBH_CF_smarr
