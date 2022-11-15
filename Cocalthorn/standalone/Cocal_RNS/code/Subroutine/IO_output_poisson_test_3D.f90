subroutine IO_output_poisson_test_3D
  use phys_constant, only : long
  use def_metric, only  : alph, psi
  use def_matter, only  : emd, rs
  use def_matter_parameter, only  : ome, ber, radi
  use def_binary_parameter, only  : dis
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  implicit none
  integer :: ir, it, ip
!
! --- Matter
  open(12,file='rnsflu_3D.las',status='unknown')
  write(12,'(5i5)') nrf, ntf, npf
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        write(12,'(1p,6e20.12)') emd(ir,it,ip), rs(it,ip)
      end do
    end do
  end do
  write(12,'(1p,6e20.12)') ome, ber, radi, dis
  close(12)
!
! --- Metric potentials.
  open(13,file='rnsgra_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ip = 0, npg
    do it = 0, ntg
      do ir = 0, nrg
        write(13,'(1p,6e20.12)')  psi(ir,it,ip), &
    &                            alph(ir,it,ip)
      end do
    end do
  end do
  write(13,'(1p,6e20.12)') ome, ber, radi, dis
  close(13)
!
! --- Coordinate grids.
  open(14,file='rnsgrids_3D.las',status='unknown')
  write(14,'(5i5)') nrg, ntg, npg
  do ir = 0, nrg
    write(14,'(1p,6e20.12)') rg(ir)
  end do
  do it = 0, ntg
    write(14,'(1p,6e20.12)') thg(it)
  end do
  do ip = 0, npg
    write(14,'(1p,6e20.12)') phig(ip)
  end do
  close(14)
!
end subroutine IO_output_poisson_test_3D
