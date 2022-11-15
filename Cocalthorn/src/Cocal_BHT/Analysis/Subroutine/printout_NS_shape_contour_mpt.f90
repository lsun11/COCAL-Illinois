subroutine printout_NS_shape_contour_mpt(impt)
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use grid_parameter, only : ntf, npf, ntfeq, ntfxy, &
  &                          npfxzp, npfxzm, npfyzp, npfyzm
  use def_matter, only : rs
  implicit none
  integer :: iseq
  integer :: it, ip, impt
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/), char_1
!
  open(20,file='bnsshape_xy_mpt'//np(impt)//'.dat',status='unknown')
  open(21,file='bnsshape_xz_mpt'//np(impt)//'.dat',status='unknown')
  open(22,file='bnsshape_yz_mpt'//np(impt)//'.dat',status='unknown')
!
  it = ntfeq
  do ip = 0, npf
    write(20,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*cosphig(ip), rs(it,ip)*sinthg(it)*sinphig(ip)
  end do
  write(20,'(/)') 
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
  write(21,'(/)') 
!
  ip = npfyzp
  do it = 0, ntf
    write(22,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*sinphig(ip), rs(it,ip)*costhg(it)
  end do
  ip = npfyzm
  do it = ntf-1, 0, -1
    write(22,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*sinphig(ip), rs(it,ip)*costhg(it)
  end do
  write(22,'(/)') 
!
end subroutine printout_NS_shape_contour_mpt
