subroutine write_last_physq
  use def_quantities
  use def_bh_parameter
  use grid_parameter, only : rgin
  implicit none

  open(15,file='last_physq.txt' ,status='unknown')
  write(15,'(1p,13e16.8)') rgin, spin_bh, ome_bh, admmass, komarmass, admmass_thr,  &
                   &       app_hor_area_bh, irredmass, bindingene,   &
                   &       angmom, angmom_thr, angmom_smarr, qua_loc_spin

  write (15,'(a30,1p,e16.8)') ' Radius of bh               = ', rgin
  write (15,'(a30,1p,e16.8)') ' Spin of                    = ', spin_bh
  write (15,'(a30,1p,e16.8)') ' Orbital angular velocity   = ', ome_bh
  write (15,'(a30,1p,e16.8)') ' ADM   mass                 = ', admmass
  write (15,'(a30,1p,e16.8)') ' Komar mass                 = ', komarmass
  write (15,'(a30,1p,e16.8)') ' ADM mass (throats)         = ', admmass_thr
  write (15,'(a30,1p,e16.8)') ' Apparent horizon area      = ', app_hor_area_bh
  write (15,'(a30,1p,e16.8)') ' Irreducible mass           = ', irredmass
  write (15,'(a30,1p,e16.8)') ' Binding energy             = ', bindingene
  write (15,'(a30,1p,e16.8)') ' Angular momentum           = ', angmom
  write (15,'(a30,1p,e16.8)') ' Angular momentum (throats) = ', angmom_thr
  write (15,'(a30,1p,e16.8)') ' Angular momentum (Smarr)   = ', angmom_smarr
  write (15,'(a30,1p,e16.8)') ' Quasi local spin BH        = ', qua_loc_spin
  write (15,*) '#'
  close(15)

end subroutine write_last_physq
