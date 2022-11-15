subroutine printout_physq_console_BBH
  use phys_constant, only  :   long, pi
  use def_matter, only : emd 
  use def_matter_parameter, only : ome, radi
  use def_quantities
  use def_bh_parameter
  use grid_parameter, only : rgin
  implicit none
!
!  write (6,'(a20,1p,e14.6)') ' Omega M and ome  = ', ome/radi*gravmass_sph
!  write (6,'(a20,1p,e14.6)') ' Radius in x-axis = ', radi
!  write (6,'(a20,1p,e14.6)') ' (p/rho)_c        = ', emd(0,0,0)
!  write (6,'(a20,1p,e14.6)') ' Rest  mass       = ', restmass

  write (6,'(a30,1p,e16.8)') ' Radius of bh               = ', rgin
  write (6,'(a30,1p,e16.8)') ' Spin of                    = ', spin_bh
  write (6,'(a30,1p,e16.8)') ' Orbital angular velocity   = ', ome_bh
  write (6,'(a30,1p,e16.8)') ' ADM   mass                 = ', admmass
  write (6,'(a30,1p,e16.8)') ' Komar mass                 = ', komarmass
  write (6,'(a30,1p,e16.8)') ' ADM mass (throats)         = ', admmass_thr
  write (6,'(a30,1p,e16.8)') ' Apparent horizon area      = ', app_hor_area_bh
  write (6,'(a30,1p,e16.8)') ' Irreducible mass           = ', irredmass
  write (6,'(a30,1p,e16.8)') ' Binding energy             = ', bindingene
  write (6,'(a30,1p,e16.8)') ' Angular momentum           = ', angmom
  write (6,'(a30,1p,e16.8)') ' Angular momentum (throats) = ', angmom_thr
  write (6,'(a30,1p,e16.8)') ' Angular momentum (Smarr)   = ', angmom_smarr
  write (6,'(a30,1p,e16.8)') ' Quasi local spin BH        = ', qua_loc_spin

!  write (6,'(a20,1p,e14.6)') ' T/|W|            = ', ToverW
end subroutine printout_physq_console_BBH
