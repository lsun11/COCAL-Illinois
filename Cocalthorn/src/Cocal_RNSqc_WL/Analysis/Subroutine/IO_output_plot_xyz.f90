subroutine IO_output_plot_xyz
  use phys_constant, only : long
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd
  use def_matter, only : emd, omef
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : ntgeq, npgxzp, npgxzm, nrg, nrf, &
  &                          npgyzp, npgyzm, ntgpolp, ntgpolm
  implicit none
  integer :: irg, itg, ipg, irf
!
! --- For surface plot
  open(12,file='plot_x.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgeq,npgxzm), &
    &                                 alph(irg,ntgeq,npgxzm), &
    &                                 bvxd(irg,ntgeq,npgxzm), &
    &                                 bvyd(irg,ntgeq,npgxzm), &
    &                                 bvzd(irg,ntgeq,npgxzm)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgeq,npgxzp), &
    &                                 alph(irg,ntgeq,npgxzp), &
    &                                 bvxd(irg,ntgeq,npgxzp), &
    &                                 bvyd(irg,ntgeq,npgxzp), &
    &                                 bvzd(irg,ntgeq,npgxzp)

  end do
  close(12)
!
!
  open(12,file='plot_y.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgeq,npgyzm), &
    &                                 alph(irg,ntgeq,npgyzm), &
    &                                 bvxd(irg,ntgeq,npgyzm), &
    &                                 bvyd(irg,ntgeq,npgyzm), &
    &                                 bvzd(irg,ntgeq,npgyzm)

  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgeq,npgyzp), &
    &                                 alph(irg,ntgeq,npgyzp), &
    &                                 bvxd(irg,ntgeq,npgyzp), &
    &                                 bvyd(irg,ntgeq,npgyzp), &
    &                                 bvzd(irg,ntgeq,npgyzp)
  end do
  close(12)
!
!
  open(12,file='plot_z.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgpolm,0), &
    &                                 alph(irg,ntgpolm,0), &
    &                                 bvxd(irg,ntgpolm,0), &
    &                                 bvyd(irg,ntgpolm,0), &
    &                                 bvzd(irg,ntgpolm,0)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgpolp,0), &
    &                                 alph(irg,ntgpolp,0), &
    &                                 bvxd(irg,ntgpolp,0), &
    &                                 bvyd(irg,ntgpolp,0), &
    &                                 bvzd(irg,ntgpolp,0)
  end do
  close(12)
!
! --- For fluid plot
  open(12,file='plot_x_matter.dat',status='unknown')
  do irf = nrf, 0, -1
    write(12,'(1p,6e20.12)') -rg(irf), emd(irf,ntgeq,npgxzm), &
    &                                 omef(irf,ntgeq,npgxzm)
  end do
  do irf = 0, nrf
    write(12,'(1p,6e20.12)')  rg(irf), emd(irf,ntgeq,npgxzp), &
    &                                 omef(irf,ntgeq,npgxzp)
!
  end do
  close(12)
!
end subroutine IO_output_plot_xyz
