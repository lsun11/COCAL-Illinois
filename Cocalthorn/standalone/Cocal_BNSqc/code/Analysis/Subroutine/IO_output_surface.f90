subroutine IO_output_surface
  use phys_constant, only : long
  use def_matter, only  : rs
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   ntf, npf
  implicit none
  integer :: it, ip
!
! --- For surface plot
  open(12,file='rnssurface.dat',status='unknown')
  do ip = 0, npf
    do it = 0, ntf
      write(12,'(1p,6e20.12)') thg(it), phig(ip), rs(it,ip)
    end do
  end do
  close(12)
!
end subroutine IO_output_surface
