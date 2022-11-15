subroutine printout_NS_shape_mpt(impt)
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use grid_parameter, only : ntf, npf, ntfeq, ntfxy, &
  &                          npfxzp, npfxzm, npfyzp, npfyzm
  use def_matter, only : rs
  implicit none
  integer :: it, ip
  integer,intent(in)  :: impt
  character(len=1) :: np(2) = (/'1', '2'/)
!
  open(20,file='rnsshape_mpt'//np(impt)//'.dat',status='unknown')
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
    write(20,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*cosphig(ip), rs(it,ip)*costhg(it)
  end do
  ip = npfxzm
  do it = ntf-1, 0, -1
    write(20,'(1p,2e14.6)') &
  & rs(it,ip)*sinthg(it)*cosphig(ip), rs(it,ip)*costhg(it)
  end do
  write(20,'(/)') 
!
  ip = npfyzp
  do it = 0, ntf
    write(20,'(1p,2e14.6)') &
  & rs(it,ip)*costhg(it), rs(it,ip)*sinthg(it)*sinphig(ip)
  end do
  ip = npfyzm
  do it = ntf-1, 0, -1
    write(20,'(1p,2e14.6)') &
  & rs(it,ip)*costhg(it), rs(it,ip)*sinthg(it)*sinphig(ip)
  end do
  close(20)
!
end subroutine printout_NS_shape_mpt
