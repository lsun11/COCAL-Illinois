subroutine initialize_AH_ks
  use phys_constant, only  : long, pi
  use grid_parameter, only : ntg, npg, ntgeq, ntgxy, &
  &                          npgxzp, npgxzm, npgyzp, npgyzm
  use def_bh_parameter, only : th_spin_bh, phi_spin_bh
  use def_kerr_schild, only : kerr_a, reh_ks
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only  : phig  
  use def_horizon, only : ahz
  integer    :: itg, ipg
  real(long) :: a2,r2,r2pa2,cosph12,cosang
!
  a2=0.0d0
  do i=1,3
    a2 = a2 + kerr_a(i)*kerr_a(i)
  end do

  r2    = reh_ks*reh_ks
  r2pa2 = r2+a2
  do ipg = 0, npg
    do itg = 0, ntg
      cosph12 = cosphig(ipg)*dcos(phi_spin_bh) + sinphig(ipg)*dsin(phi_spin_bh)
      cosang  = sinthg(itg)*dsin(th_spin_bh)*cosph12 + costhg(itg)*dcos(th_spin_bh)

!     AH radial r in spherical code coordinates       
      ahz(itg,ipg) = dsqrt(r2*r2pa2/(r2+a2*cosang**2))        
    end do
  end do
!
  open(20,file='ah_ks_xy.dat',status='unknown')
  open(21,file='ah_ks_xz.dat',status='unknown')
  open(22,file='ah_ks_yz.dat',status='unknown')
  open(24,file='ah_ks.dat',   status='unknown')

  itg = ntgeq
  do ipg = 0, npg
    write(20,'(1p,2e14.6)') &
  & ahz(itg,ipg)*sinthg(itg)*cosphig(ipg), ahz(itg,ipg)*sinthg(itg)*sinphig(ipg)
  end do
  write(20,'(1x)')
!
  ipg = npgxzp
  do itg = 0, ntg
    write(21,'(1p,2e14.6)') &
  & ahz(itg,ipg)*sinthg(itg)*cosphig(ipg), ahz(itg,ipg)*costhg(itg)
  end do
  ipg = npgxzm
  do itg = ntg-1, 0, -1
    write(21,'(1p,2e14.6)') &
  & ahz(itg,ipg)*sinthg(itg)*cosphig(ipg), ahz(itg,ipg)*costhg(itg)
  end do
  write(21,'(1x)')
!
  ipg = npgyzp
  do itg = 0, ntg
    write(22,'(1p,2e14.6)') &
  & ahz(itg,ipg)*sinthg(itg)*sinphig(ipg), ahz(itg,ipg)*costhg(itg)
  end do
  ipg = npgyzm
  do itg = ntg-1, 0, -1
    write(22,'(1p,2e14.6)') &
  & ahz(itg,ipg)*sinthg(itg)*sinphig(ipg), ahz(itg,ipg)*costhg(itg)
  end do
  write(22,'(1x)')
!
  do ipg = 0, npg
    do itg = 0, ntg
      write(24,'(1p,6e20.12)') phig(ipg), thg(itg)-0.5d0*pi, ahz(itg,ipg)
    end do
    write(24,'(1x)')
  end do
!
  close(20)
  close(21)
  close(22)
  close(24)
!
end subroutine initialize_AH_ks
