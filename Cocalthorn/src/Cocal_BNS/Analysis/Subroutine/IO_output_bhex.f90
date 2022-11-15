subroutine IO_output_bhex
  use phys_constant, only : long
  use def_metric, only : psi, bvxd
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : ntgeq, npgxzp, npgxzm, nrg, npgyzp, npgyzm, ntgpolp, ntgpolm
  implicit none
  integer :: irg, itg, ipg
!
! --- For surface plot
  open(12,file='plot_x.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgeq,npgxzm) &
    &                                , bvxd(irg,ntgeq,npgxzm)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgeq,npgxzp) &
    &                                , bvxd(irg,ntgeq,npgxzp)
  end do
  close(12)
!
!
  open(12,file='plot_y.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgeq,npgyzm) &
    &                                , bvxd(irg,ntgeq,npgyzm)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgeq,npgyzp) &
    &                                , bvxd(irg,ntgeq,npgyzp)
  end do
  close(12)
!
!
  open(12,file='plot_z.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgpolm,0) &
    &                                , bvxd(irg,ntgpolm,0)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgpolp,0) &
    &                                , bvxd(irg,ntgpolp,0)
  end do
  close(12)
!
end subroutine IO_output_bhex
