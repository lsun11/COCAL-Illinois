subroutine calc_ang_mom_peos_irrot
!
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrf, ntf, npf
  use def_matter_parameter, only : radi
  use def_quantities, only : angmom
  use make_array_3d
  use interface_source_ang_mom_peos_irrot
  use interface_vol_int_fluid
  implicit none
  real(long)          :: volf
  real(long), pointer :: souf(:,:,:)
  integer :: ir, it, ip
!
  call alloc_array3d(souf,0, nrf,0, ntf,0, npf)
!
  call source_ang_mom_peos_irrot(souf)
  call vol_int_fluid(souf,volf)
  angmom = radi**4*volf
!
  write (6,'(a20,1p,e14.6)') ' Angular momentum = ', angmom
!
  deallocate(souf)
!
end subroutine calc_ang_mom_peos_irrot
