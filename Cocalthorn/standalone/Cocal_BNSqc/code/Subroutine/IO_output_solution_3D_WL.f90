subroutine IO_output_solution_3D_WL
  use phys_constant, only : long
  use def_metric_hij, only  : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
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
        write(13,'(1p,6e20.12)') hxxd(irg,itg,ipg), &
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
end subroutine IO_output_solution_3D_WL
