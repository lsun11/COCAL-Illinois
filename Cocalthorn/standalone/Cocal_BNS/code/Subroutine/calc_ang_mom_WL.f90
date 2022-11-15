subroutine calc_ang_mom_WL
!
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter_parameter, only : radi
  use def_quantities, only : angmom
  use make_array_3d
  use interface_source_ang_mom_WL
  use interface_vol_int_grav
  use interface_vol_int_fluid
  implicit none
  real(long)          :: volf, volg, fac8pi
  real(long), pointer :: soug(:,:,:), souf(:,:,:)
!
  call alloc_array3d(soug, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(souf,0, nrf,0, ntf,0, npf)
!
  call source_ang_mom_WL(soug,souf)
  call vol_int_grav(soug,volg)
  call vol_int_fluid(souf,volf)
!
  fac8pi = 0.125d0/pi
  angmom = radi**4*volf + fac8pi*radi**2*volg
!
!  write (6,'(a20,1p,e14.6)') ' Angular momentum = ', angmom
!
  deallocate(soug)
  deallocate(souf)
!
end subroutine calc_ang_mom_WL
