subroutine IO_input_initial_1D
  use phys_constant, only : long
  use grid_parameter
  use coordinate_grav_r, only : rg
  use def_metric
  use def_matter
  use def_matter_parameter, only  : ome, ber, radi
  implicit none
  real(long)    ::     tmp, tmp1, dummy
  integer       ::     nrftmp, ir, nrgtmp
!
! --- Matter variables.
!
  rs = 1.0d0
  emdg = 0.0d0
  ome = 1.0d-2
  omef(0:nrf, 0:ntf, 0:npf) = ome
  omeg(0:nrf, 0:ntf, 0:npf) = ome
!
  open(2,file='rnsflu_1D.ini',status='old')
!
  read(2,'(5i5)') nrftmp
  do ir = 0, nrftmp
    read(2,'(1p,6e20.12)') dummy, tmp
    emdg(ir, 0:ntg, 0:npg) = tmp
  end do
  emd(0:nrf, 0:ntf, 0:npf) = emdg(0:nrf, 0:ntf, 0:npf)
  read(2,'(1p,6e20.12)') ber, radi
!
  close(2)
!
! --- Metric potentials.
!
  open(3,file='rnsgra_1D.ini',status='old')
!
  read(3,'(5i5)') nrgtmp
  do ir = 0, nrgtmp
    read(3,'(1p,6e20.12)') dummy, tmp, tmp1
    psi(ir, 0:ntg, 0:npg) = tmp
    alph(ir, 0:ntg, 0:npg) = tmp1
    alps = tmp * tmp1
  end do
!
  close(3)
  bvxd = 0.0d0
  bvyd = 0.0d0
  bvzd = 0.0d0
  tfkij = 0.0d0
  tfkijkij = 0.0d0
end subroutine IO_input_initial_1D
