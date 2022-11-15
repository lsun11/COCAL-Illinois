subroutine IO_output_plot_xyz_BBH_CF_mpt(impt)
  use phys_constant, only : long
  use def_metric
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : ntgeq, npgxzp, npgxzm, nrg, &
  &                           npgyzp, npgyzm, ntgpolp, ntgpolm
  use grid_points_binary_excision, only : irg_exin, irg_exout
  use grid_parameter_binary_excision, only :ex_radius 
  use def_binary_parameter, only : sepa, dis
  implicit none
  integer :: irg, itg, ipg, impt, npg_l, npg_r, count
  real(long) :: work_shift, pari
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/), char_1
!
  if (impt.eq.1) then
    npg_l = npgxzm
    npg_r = npgxzp
    work_shift = - dis
    pari = 1.0
  else if(impt.eq.2) then
    npg_l = npgxzp
    npg_r = npgxzm
    work_shift = dis
    pari = -1.0
  else if(impt.eq.3) then
    npg_l = npgxzm
    npg_r = npgxzp
    work_shift = 0.0d0
    pari = 1.0
  end if
!
  open(12,file='plot_x_mpt'//np(impt)//'.dat',status='unknown')
  count = 0
  do irg = nrg, 0, -1
    char_1 = ' '
    if (impt.ne.3) then
      if (irg.gt.irg_exin(ntgeq,npg_l).and.irg.lt.irg_exout(ntgeq,npg_l)) then
        char_1 = '#'
        count = count + 1
      end if
    end if
    write(12,'(a1,1p,10e20.12)') char_1,  -rg(irg)+work_shift   &
    &                                  ,  psi(irg,ntgeq,npg_l) &
    &                                  , alph(irg,ntgeq,npg_l) &
    &                                  , bvxd(irg,ntgeq,npg_l)*pari &
    &                                  , bvyd(irg,ntgeq,npg_l)*pari &
    &                                  , bvzd(irg,ntgeq,npg_l)
    if (count.eq.1) write(12,'(/)')
  end do
  write(12,'(/)')
  count = 0
  do irg = 0, nrg
    char_1 = ' '
    if (impt.ne.3) then
      if (irg.gt.irg_exin(ntgeq,npg_r).and.irg.lt.irg_exout(ntgeq,npg_r)) then
        char_1 = '#'
        count = count + 1
      end if
    end if
    write(12,'(a1,1p,10e20.12)') char_1,   rg(irg)+work_shift   &
      &                                ,  psi(irg,ntgeq,npg_r) &
      &                                , alph(irg,ntgeq,npg_r) &
      &                                , bvxd(irg,ntgeq,npg_r)*pari &
      &                                , bvyd(irg,ntgeq,npg_r)*pari &
      &                                , bvzd(irg,ntgeq,npg_r)
    if (count.eq.1) write(12,'(/)')
  end do
  close(12)
!
!
  open(12,file='plot_y_mpt'//np(impt)//'.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg),  psi(irg,ntgeq,npgyzm)  &
    &                                , alph(irg,ntgeq,npgyzm)  &
    &                                , bvxd(irg,ntgeq,npgyzm)  &
    &                                , bvyd(irg,ntgeq,npgyzm)  &
    &                                , bvzd(irg,ntgeq,npgyzm)
  end do
  write(12,*) ' '
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg),  psi(irg,ntgeq,npgyzp)  &
    &                                , alph(irg,ntgeq,npgyzp)  &
    &                                , bvxd(irg,ntgeq,npgyzp)  &
    &                                , bvyd(irg,ntgeq,npgyzp)  &
    &                                , bvzd(irg,ntgeq,npgyzp)
  end do
  close(12)
!
!
  open(12,file='plot_z_mpt'//np(impt)//'.dat',status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,6e20.12)') -rg(irg),  psi(irg,ntgpolm,0)  &
    &                                , alph(irg,ntgpolm,0)  &
    &                                , bvxd(irg,ntgpolm,0)  &
    &                                , bvyd(irg,ntgpolm,0)  &
    &                                , bvzd(irg,ntgpolm,0)
  end do
  write(12,*) ' '
  do irg = 0, nrg
    write(12,'(1p,6e20.12)')  rg(irg),  psi(irg,ntgpolp,0) &
    &                                , alph(irg,ntgpolp,0) &
    &                                , bvxd(irg,ntgpolp,0) &
    &                                , bvyd(irg,ntgpolp,0) &
    &                                , bvzd(irg,ntgpolp,0)
  end do
  close(12)
!
end subroutine IO_output_plot_xyz_BBH_CF_mpt
