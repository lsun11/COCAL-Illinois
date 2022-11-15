subroutine IO_printout_grid_data_noex
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg, drg
  use coordinate_grav_theta
  use coordinate_grav_phi
  implicit none
  integer :: irg
  open(1,file='grid_data.dat',status='unknown')
  write(1,'(a4,i3,a10,1p,e20.12)') 'ntg=', ntg, '     dthg=', dthg
  write(1,'(a4,i3,a10,1p,e20.12)') 'npg=', npg, '     dphg=', dphig
!
  write(1,'(a37,i3,a2,1p,e20.12)')  '..................................rg(',0,')=', rg(0)
  do irg=1, nrg
    write(1,'(a4,i3,a2,1p,e20.12,a5,a3,i3,a2,1p,e20.12)') 'drg(',irg,')=',drg(irg),'     ','rg(',irg,')=',rg(irg)
  end do
  close(1)
end subroutine IO_printout_grid_data_noex

