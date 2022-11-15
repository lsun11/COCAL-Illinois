subroutine IO_output_solution_3D_CF_BH
  use phys_constant, only : long
  use def_metric, only  : alph, psi, bvxd, bvyd, bvzd
  use def_matter, only  : emd, rs
  use def_bh_parameter, only : ome_bh, spin_bh, alph_bh
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   nrg, ntg, npg
  implicit none
  integer :: ir, it, ip
!
! --- Metric potentials.
  open(13,file='bbhgra_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ip = 0, npg
    do it = 0, ntg
      do ir = 0, nrg
        write(13,'(1p,6e20.12)')  psi(ir,it,ip), &
    &                            alph(ir,it,ip), &
    &                            bvxd(ir,it,ip), &
    &                            bvyd(ir,it,ip), &
    &                            bvzd(ir,it,ip)
      end do
    end do
  end do
  write(13,'(1p,6e20.12)') ome_bh, spin_bh, alph_bh
  close(13)
!
! --- Coordinate grids.
  open(14,file='bbhgrids_3D.las',status='unknown')
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
end subroutine IO_output_solution_3D_CF_BH
