subroutine read_omega
  use phys_constant, only : long, pi
  use def_bh_parameter, only : ome_bh
  implicit none
  open(1,file='omega_last.txt',status='old')
  read(1,'(1p,e20.12)') ome_bh
  close(1)
end subroutine read_omega
