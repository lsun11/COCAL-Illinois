subroutine IO_output_solution_3D_WL_BHT
  use phys_constant, only : long
  use def_metric_hij
!  use def_metric_excurve_grid, only : trk_grid
!  use def_metric_dihiju, only : dihixu_grid, dihiyu_grid, dihizu_grid
  use def_CTT_decomposition, only : wvxd, wvyd, wvzd, sigt
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   nrg, ntg, npg
  implicit none
  integer :: irg, itg, ipg
!
! --- Metric potentials.
  open(13,file='rnsgra_hij_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        write(13,'(1p,6e23.15)') hxxd(irg,itg,ipg), &
        &                        hxyd(irg,itg,ipg), &
        &                        hxzd(irg,itg,ipg), &
        &                        hyyd(irg,itg,ipg), &
        &                        hyzd(irg,itg,ipg), &
        &                        hzzd(irg,itg,ipg)
      end do
    end do
  end do
  close(13)
!
! --- Gauge potentials.
  open(14,file='rnsgra_gauge_3D.las',status='unknown')
  write(14,'(5i5)') nrg, ntg, npg
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        write(14,'(1p,10e23.15)')   gaugex(irg,itg,ipg), &
        &                           gaugey(irg,itg,ipg), &
        &                           gaugez(irg,itg,ipg), &
        &                             wvxd(irg,itg,ipg), &
        &                             wvyd(irg,itg,ipg), &
        &                             wvzd(irg,itg,ipg), &
        &                             sigt(irg,itg,ipg)
      end do
    end do
  end do
  close(14)

end subroutine IO_output_solution_3D_WL_BHT
