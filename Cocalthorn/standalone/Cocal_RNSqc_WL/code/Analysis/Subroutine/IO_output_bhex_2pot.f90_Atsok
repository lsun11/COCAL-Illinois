subroutine IO_output_bhex_2pot
  use phys_constant, only : long
  use def_metric, only : psi, alph, bvxd, bvyd
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : ntgeq, npgxzp, npgxzm, nrg, npgyzp, npgyzm, ntgpolp, ntgpolm
  use grid_points_binary_excision, only : rb
  use grid_parameter_binary_excision, only :ex_radius 
  implicit none
  integer :: irg, itg, ipg
!
! --- For surface plot
  open(12,file='plot_x.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgeq,npgxzm)  &
    &                                , bvxd(irg,ntgeq,npgxzm) &
    &                                , alph(irg,ntgeq,npgxzm) &
    &                                , bvyd(irg,ntgeq,npgxzm)
  end do
  write(12,*) ' '
  do irg = 0, nrg
!    if (rb(irg,ntgeq,npgxzp).ge.rg(0)) then
    if (rb(irg,ntgeq,npgxzp).ge.ex_radius) then
      write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgeq,npgxzp)  &
      &                                , bvxd(irg,ntgeq,npgxzp) &
      &                                , alph(irg,ntgeq,npgxzp) &
      &                                , bvyd(irg,ntgeq,npgxzp)
    else
      write(12,'(a1,1p,6e20.12)')  '#', rg(irg), psi(irg,ntgeq,npgxzp)  &
      &                                        , bvxd(irg,ntgeq,npgxzp) &
      &                                        , alph(irg,ntgeq,npgxzp) &
      &                                        , bvyd(irg,ntgeq,npgxzp)
    end if 
  end do
  close(12)
!
!
  open(12,file='plot_y.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgeq,npgyzm)   &
    &                                , bvxd(irg,ntgeq,npgyzm)  &
    &                                , alph(irg,ntgeq,npgyzm)  &
    &                                , bvyd(irg,ntgeq,npgyzm)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgeq,npgyzp)   &
    &                                , bvxd(irg,ntgeq,npgyzp)  &
    &                                , alph(irg,ntgeq,npgyzp)  &
    &                                , bvyd(irg,ntgeq,npgyzp)
  end do
  close(12)
!
!
  open(12,file='plot_z.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg), psi(irg,ntgpolm,0)   &
    &                                , bvxd(irg,ntgpolm,0)  &
    &                                , alph(irg,ntgpolm,0)  &
    &                                , bvyd(irg,ntgpolm,0)
  end do
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg), psi(irg,ntgpolp,0)  &
    &                                , bvxd(irg,ntgpolp,0) &
    &                                , alph(irg,ntgpolp,0) &
    &                                , bvyd(irg,ntgpolp,0)
  end do
  close(12)
!
end subroutine IO_output_bhex_2pot
