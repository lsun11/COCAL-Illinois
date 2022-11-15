subroutine printout_debug
  use coordinate_grav_r, only : rg
  use grid_parameter, only : nrg, ntg, npg, ntf, npf, ntfeq, ntfxy, &
  &                          npfxzp, npfxzm, npfyzp, npfyzm
  use def_matter, only : rs
  use def_metric, only : psi, alph, bvyd
  implicit none
  integer :: it, ip, irg, itg, ipg
!
!!  open(82,file='plot_metric.dat',status='unknown')
!!    itg = ntg/2
!!    ipg = 0
!!    do irg = 0, nrg
!!      write(82,'(1p,10e14.6)') rg(irg), psi(irg,itg,ipg), & 
!!      &             alph(irg,itg,ipg), -bvyd(irg,itg,ipg)
!!!!      &        souvec(irg,itg,ipg,1), &
!!!!      &        souvec(irg,itg,ipg,2), &
!!!!      &        souvec(irg,itg,ipg,3), soutest(irg,itg,ipg)
!!    end do
!!  close(82)
!
  write(6,*)       ' Coordinate radius in xz plane '
  do it = 0, ntf
    write(6,'(i5,1p,e14.6)') it, rs(it,npfxzp)
  end do
  write(6,*)       ' Coordinate radius in yz plane '
  do it = 0, ntf
    write(6,'(i5,1p,e14.6)') it, rs(it,npfyzp)
  end do
  write(6,*)       ' Coordinate radius in xy plane '
  do ip = 0, npf
    write(6,'(i5,1p,e14.6)') ip, rs(ntfeq,ip)
  end do
end subroutine printout_debug
