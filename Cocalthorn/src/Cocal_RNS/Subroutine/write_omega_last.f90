subroutine write_omega_last
  use def_bh_parameter, only : ome_bh
  use def_quantities
  use grid_parameter, only : rgin
  implicit none

  open(15,file='omega_last.txt' ,status='unknown')
  write(15,'(1p,e10.4,1p,10e16.8)') rgin, ome_bh, admmass, komarmass, admmass_thr,  &
                   &                 app_hor_area_bh, irredmass, bindingene,   &
                   &                 angmom, angmom_thr, angmom_smarr
  close(15)

end subroutine write_omega_last
