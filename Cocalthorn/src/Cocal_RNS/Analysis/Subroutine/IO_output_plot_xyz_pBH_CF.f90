subroutine IO_output_plot_xyz_pBH_CF
  use phys_constant, only : long
  use def_metric, only : psi, alph
  use def_metric_pBH, only : psi_trpBH, alph_trpBH, wme, index_wme
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : ntgeq, npgxzp, npgxzm, nrg, &
  &                          npgyzp, npgyzm, ntgpolp, ntgpolm
  implicit none
  integer :: irg, itg, ipg
!
! --- For surface plot
  open(12,file='plot_x.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,20e20.12)')-rg(irg), psi(irg,ntgeq,npgxzm), psi_trpBH(irg), &
    &                                 alph(irg,ntgeq,npgxzm),alph_trpBH(irg), &
    &                                  wme(irg,ntgeq,npgxzm), psi_trpBH(irg) &
    &                                                         **(-index_wme)
  end do
  do irg = 0, nrg
    write(12,'(1p,20e20.12)') rg(irg), psi(irg,ntgeq,npgxzp), psi_trpBH(irg), &
    &                                 alph(irg,ntgeq,npgxzp),alph_trpBH(irg), &
    &                                  wme(irg,ntgeq,npgxzp), psi_trpBH(irg) &
    &                                                         **(-index_wme)
  end do
  close(12)
!
!
  open(12,file='plot_y.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,20e20.12)')-rg(irg), psi(irg,ntgeq,npgyzm), psi_trpBH(irg), &
    &                                 alph(irg,ntgeq,npgyzm),alph_trpBH(irg), &
    &                                  wme(irg,ntgeq,npgyzm), psi_trpBH(irg) &
    &                                                         **(-index_wme)
  end do
  do irg = 0, nrg
    write(12,'(1p,20e20.12)') rg(irg), psi(irg,ntgeq,npgyzp), psi_trpBH(irg), &
    &                                 alph(irg,ntgeq,npgyzp),alph_trpBH(irg), &
    &                                  wme(irg,ntgeq,npgyzp), psi_trpBH(irg) &
    &                                                         **(-index_wme)
  end do
  close(12)
!
!
  open(12,file='plot_z.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,20e20.12)')-rg(irg), psi(irg,ntgpolm,0), psi_trpBH(irg), &
    &                                 alph(irg,ntgpolm,0),alph_trpBH(irg), &
    &                                  wme(irg,ntgpolm,0), psi_trpBH(irg) &
    &                                                      **(-index_wme)
  end do
  do irg = 0, nrg
    write(12,'(1p,20e20.12)') rg(irg), psi(irg,ntgpolp,0), psi_trpBH(irg), &
    &                                 alph(irg,ntgpolp,0),alph_trpBH(irg), &
    &                                  wme(irg,ntgpolp,0), psi_trpBH(irg) &
    &                                                      **(-index_wme)
  end do
  close(12)
!
end subroutine IO_output_plot_xyz_pBH_CF
