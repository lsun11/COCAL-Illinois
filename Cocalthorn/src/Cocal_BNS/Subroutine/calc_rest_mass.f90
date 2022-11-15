subroutine calc_rest_mass
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrf, ntf, npf
  use def_matter_parameter, only : radi
  use make_array_3d
  use def_quantities, only : restmass
  use interface_source_rest_mass
  use interface_vol_int_fluid
  implicit none
  real(long)          :: fac2pi, fac4pi
  real(long)          :: volf
  real(long), pointer :: souf(:,:,:)
  integer :: ir, it, ip
!
  call alloc_array3d(souf, 0, nrf, 0, ntf, 0, npf)
!      
  call source_rest_mass(souf)
  call vol_int_fluid(souf,volf)
!
  restmass = radi**3*volf
!
!  write (6,'(a20,1p,e14.6)') ' Rest  mass =       ',  restmass
!
  deallocate(souf)
end subroutine calc_rest_mass
