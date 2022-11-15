subroutine IO_output_solution_1D
  use phys_constant, only : long
  use def_metric, only  : alph, psi
  use def_matter
  use def_matter_parameter, only : ber, radi
  use coordinate_grav_r, only : rg
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf
  implicit none
  integer :: ir, it
!
! For spherical
  open(12,file='rnsflu_1D.las',status='unknown')
  write(12,'(5i5)') nrf
  do ir = 0, nrf
    write(12,'(1p,6e20.12)') rg(ir), emdg(ir, ntg, 0)
  end do
  write(12,'(1p,6e20.12)') ber, radi
  close(12)
  do it = 0, ntf
  write(6,*)rs(it,0), emd(nrf,it,0)
  end do
!
! --- Metric potentials.
  open(13,file='rnsgra_1D.las',status='unknown')
  write(13,'(5i5)') nrg
  do ir = 0, nrg
    write(13,'(1p,6e20.12)') rg(ir), psi(ir, ntg, 0), &
    &                               alph(ir, ntg, 0), &
    &                               emdg(ir, ntg, 0)
  end do
  write(13,'(1p,6e20.12)') ber, radi
  close(13)
!
end subroutine IO_output_solution_1D
