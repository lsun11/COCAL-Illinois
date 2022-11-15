subroutine printout_axis_spin_BNS_mpt(impt,filename)
  use phys_constant, only : long
  use def_metric
  use def_matter
  use def_velocity_potential
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
  character(30) :: filename
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
  open(12,file=filename,status='unknown')
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
    &                                  ,  emd(irg,ntgeq,npg_l)  &
    &                                  ,  vep(irg,ntgeq,npg_l)  &
    &                                  ,  vepxf(irg,ntgeq,npg_l)*pari  &
    &                                  ,  vepyf(irg,ntgeq,npg_l)*pari
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
      &                                ,  emd(irg,ntgeq,npg_r)   &
      &                                ,  vep(irg,ntgeq,npg_r)  &
      &                                ,  vepxf(irg,ntgeq,npg_r)*pari  &
      &                                ,  vepyf(irg,ntgeq,npg_r)*pari
    if (count.eq.1) write(12,'(/)')
  end do
  close(12)
!
end subroutine printout_axis_spin_BNS_mpt
