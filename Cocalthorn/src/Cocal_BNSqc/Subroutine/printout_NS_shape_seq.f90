subroutine printout_NS_shape_seq(iseq)
  use phys_constant, only : long, pi
  use coordinate_grav_theta, only  : thg
  use coordinate_grav_phi, only  : phig
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use grid_parameter, only : ntf, npf, ntfeq, ntfxy, &
  &                          npfxzp, npfxzm, npfyzp, npfyzm
  use def_matter, only : rs
  implicit none
  integer :: iseq
  integer :: it, ip
!
  if (iseq.eq.1) then
    open(20,file='rnsshape_seq_xy.dat',status='unknown')
    open(21,file='rnsshape_seq_xz.dat',status='unknown')
    open(22,file='rnsshape_seq_zy.dat',status='unknown')
    open(23,file='rnsshape_seq_yz.dat',status='unknown')
    open(24,file='rnsshape_NS_gnuplot.dat',status='unknown')
  else
    open(20,file='rnsshape_seq_xy.dat',status='old', position="append")
    open(21,file='rnsshape_seq_xz.dat',status='old', position="append")
    open(22,file='rnsshape_seq_zy.dat',status='old', position="append")
    open(23,file='rnsshape_seq_yz.dat',status='old', position="append")
    open(24,file='rnsshape_NS_gnuplot.dat',status='old', position="append")
  end if
!
  it = ntfeq
  do ip = 0, npf
    write(20,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*cosphig(ip), rs(it,ip)*sinthg(it)*sinphig(ip)
  end do
  write(20,'(1x)') 
!
  ip = npfxzp
  do it = 0, ntf
    write(21,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*cosphig(ip), rs(it,ip)*costhg(it)
  end do
  ip = npfxzm
  do it = ntf-1, 0, -1
    write(21,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*cosphig(ip), rs(it,ip)*costhg(it)
  end do
  write(21,'(1x)') 
!
  ip = npfyzp
  do it = 0, ntf
    write(22,'(1p,2e14.6)') &
  & rs(it,ip)*costhg(it), rs(it,ip)*sinthg(it)*sinphig(ip)
  end do
  ip = npfyzm
  do it = ntf-1, 0, -1
    write(22,'(1p,2e14.6)') &
  & rs(it,ip)*costhg(it), rs(it,ip)*sinthg(it)*sinphig(ip)
  end do
  write(22,'(1x)') 
!
  ip = npfyzp
  do it = 0, ntf
    write(23,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*sinphig(ip), rs(it,ip)*costhg(it)
  end do
  ip = npfyzm
  do it = ntf-1, 0, -1
    write(23,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*sinphig(ip), rs(it,ip)*costhg(it)
  end do
  write(23,'(1x)')
!
  do ip = 0, npf
    do it = 0, ntf
      write(24,'(1p,6e20.12)') phig(ip), thg(it)-0.5d0*pi, rs(it,ip)
    end do
    write(24,'(1x)')
  end do

  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
!
end subroutine printout_NS_shape_seq
