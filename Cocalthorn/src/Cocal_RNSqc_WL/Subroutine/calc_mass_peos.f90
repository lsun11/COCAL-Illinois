subroutine calc_mass_peos
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg, nrf, ntf, npf
  use def_matter_parameter, only : radi
  use make_array_3d
  use def_quantities, only : admmass, komarmass, komarmass_nc
  use interface_source_adm_mass_peos
  use interface_source_komar_mass_peos
  use interface_source_komar_mass_compact_peos
  use interface_vol_int_grav
  use interface_vol_int_fluid
  implicit none
  real(long)       ::     fac2pi, fac4pi
  real(long)       ::     volg, volf
  real(long),pointer :: soug(:,:,:), souf(:,:,:)
!
!
  call alloc_array3d(soug, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(souf, 0, nrf, 0, ntf, 0, npf)
!      
  call source_adm_mass_peos(soug,souf)
  call vol_int_grav(soug,volg)
  call vol_int_fluid(souf,volf)
!
  fac2pi = 0.5d0/pi
  admmass = fac2pi*(radi*volg + radi**3*volf)
!
  call source_komar_mass_peos(soug,souf)
  call vol_int_grav(soug,volg) 
  call vol_int_fluid(souf,volf)
!
  fac4pi = 0.25d0/pi
  komarmass_nc = fac4pi*(radi*volg + radi**3*volf)
!
  call source_komar_mass_compact_peos(souf)
  call vol_int_fluid(souf,volf)
!
  komarmass = radi**3*volf
!
  write (6,'(a20,1p,e14.6)') ' ADM   mass =       ', admmass
  write (6,'(a20,1p,e14.6)') ' Komar mass compact=', komarmass
  write (6,'(a20,1p,e14.6)') ' Komar mass noncomp=', komarmass_nc
!
  deallocate(soug)
  deallocate(souf)
end subroutine calc_mass_peos
