subroutine IO_output_solution_3D_BHT
  use phys_constant, only : long
  use def_metric, only  : alph, psi, bvxd, bvyd, bvzd
  use def_matter, only  : emdg, omeg
  use def_matter_parameter, only  : ome, ber, radi, emdc
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  implicit none
  integer :: ir, it, ip
!
! --- Matter
  open(12,file='rnsflu_3D.las',status='unknown')
  write(12,'(5i5)') nrg, ntg, npg
  do ip = 0, npg
    do it = 0, ntg
      do ir = 0, nrg
        write(12,'(1p,6e23.15)') emdg(ir,it,ip), omeg(ir,it,ip)
      end do
    end do
  end do
  write(12,'(1p,6e23.15)') ome, ber, radi
  write(12,'(1p,6e23.15)') emdc
  close(12)
!
! --- Metric potentials.
  open(13,file='rnsgra_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ip = 0, npg
    do it = 0, ntg
      do ir = 0, nrg
        write(13,'(1p,6e23.15)')  psi(ir,it,ip), &
    &                            alph(ir,it,ip), &
    &                            bvxd(ir,it,ip), &
    &                            bvyd(ir,it,ip), &
    &                            bvzd(ir,it,ip)
      end do
    end do
  end do
  write(13,'(1p,6e23.15)') ome, ber, radi
  close(13)
!
! --- Coordinate grids.
  open(14,file='rnsgrids_3D.las',status='unknown')
  write(14,'(5i5)') nrg, ntg, npg
  do ir = 0, nrg
    write(14,'(1p,6e23.15)') rg(ir)
  end do
  do it = 0, ntg
    write(14,'(1p,6e23.15)') thg(it)
  end do
  do ip = 0, npg
    write(14,'(1p,6e23.15)') phig(ip)
  end do
  close(14)
!
end subroutine IO_output_solution_3D_BHT
