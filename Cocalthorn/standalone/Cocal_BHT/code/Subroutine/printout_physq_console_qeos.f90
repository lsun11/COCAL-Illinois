subroutine printout_physq_console_qeos
  use def_matter, only : rhof
  use def_matter_parameter, only : ome, radi
  use def_quantities
  implicit none
!
  write (6,'(a20,1p,e14.6)') ' Omega M and ome  = ', ome/radi*gravmass_sph
  write (6,'(a20,1p,e14.6)') ' Radius in x-axis = ', radi
  write (6,'(a20,1p,e14.6)') ' rho_c            = ', rhof(0,0,0)
  write (6,'(a20,1p,e14.6)') ' Rest  mass       = ', restmass
  write (6,'(a20,1p,e14.6)') ' ADM   mass       = ', admmass
  write (6,'(a20,1p,e14.6)') ' Komar mass       = ', komarmass
  write (6,'(a20,1p,e14.6)') ' Angular momentum = ', angmom
  write (6,'(a20,1p,e14.6)') ' Charge           = ', charge
  write (6,'(a20,1p,e14.6)') ' T/|W|            = ', ToverW
end subroutine printout_physq_console_qeos
