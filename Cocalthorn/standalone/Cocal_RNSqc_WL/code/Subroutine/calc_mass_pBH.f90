subroutine calc_mass_pBH
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg
  use make_array_3d
  use def_quantities, only : admmass, komarmass, komarmass_nc
  use def_metric_pBH, only : index_wme
  use interface_sourceterm_HaC_CF_pBH
  use interface_sourceterm_trG_CF_pBH
  use interface_vol_int_grav
  implicit none
  real(long)       ::     fac2pin, fac4pi
  real(long)       ::     volg, volf
  real(long),pointer :: soug(:,:,:)
!
!
  call alloc_array3d(soug, 0, nrg, 0, ntg, 0, npg)
!      
  call sourceterm_HaC_CF_pBH(soug)
  call vol_int_grav(soug,volg)
!
  fac2pin = 0.5d0/(pi*dble(index_wme))
  admmass = fac2pin*volg
!
  call sourceterm_trG_CF_pBH(soug)
  call vol_int_grav(soug,volg)
  fac4pi = 0.25d0/pi
  komarmass = sqrt(2.0d0)*fac4pi*volg
!
!  write (6,'(a20,1p,e14.6)') ' ADM   mass =       ', admmass
!  write (6,'(a20,1p,e14.6)') ' Komar mass =       ', komarmass
!
  deallocate(soug)
end subroutine calc_mass_pBH
