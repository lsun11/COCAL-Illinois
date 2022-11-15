subroutine IO_output_surface_BNS_mpt(impt)
  use phys_constant, only : long
  use def_matter, only  : rs
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   ntf, npf
  implicit none
  integer :: it, ip, impt
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/), char_1
!
! --- For surface plot
  open(12,file='bnssurface_mpt'//np(impt)//'.dat',status='unknown')
  do ip = 0, npf
    do it = 0, ntf
      write(12,'(1p,6e20.12)') thg(it), phig(ip), rs(it,ip)
    end do
  end do
  close(12)
!
end subroutine IO_output_surface_BNS_mpt
