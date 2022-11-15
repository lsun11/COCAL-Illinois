subroutine printout_axis_CF(filename)
  use phys_constant, only : long
  use def_metric
  use def_matter
  use coordinate_grav_r, only : rg
  use grid_parameter, only  : ntgeq, npgxzp, npgxzm, nrg, &
  &                           npgyzp, npgyzm, ntgpolp, ntgpolm
  implicit none
  integer :: irg, itg, ipg
  character(30) :: filename
!
  open(12,file=filename,status='unknown')
  do irg = nrg, 0, -1
    write(12,'(1p,10e20.12)')     -rg(irg)   &
    &                           ,  psi(irg,ntgeq,npgxzm) &
    &                           , alph(irg,ntgeq,npgxzm) &
    &                           , bvxd(irg,ntgeq,npgxzm) &
    &                           , bvyd(irg,ntgeq,npgxzm) &
    &                           , emdg(irg,ntgeq,npgxzm)  
  end do
  do irg = 0, nrg
    write(12,'(1p,10e20.12)')      rg(irg)   &
      &                         ,  psi(irg,ntgeq,npgxzp) &
      &                         , alph(irg,ntgeq,npgxzp) &
      &                         , bvxd(irg,ntgeq,npgxzp) &
      &                         , bvyd(irg,ntgeq,npgxzp) &
      &                         , emdg(irg,ntgeq,npgxzp)  
  end do
  close(12)
!
end subroutine printout_axis_CF
