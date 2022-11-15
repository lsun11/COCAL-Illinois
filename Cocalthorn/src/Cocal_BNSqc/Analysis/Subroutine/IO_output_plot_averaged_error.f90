subroutine IO_output_plot_averaged_error
  use phys_constant, only : long
  use def_metric, only : psi, alph, bvxd, bvyd
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : nrg, ntg, npg
  use grid_points_binary_excision, only : irg_exin, irg_exout
  use def_binary_parameter, only : sepa, dis
  implicit none
  integer :: irg, itg, ipg, count
  real(long) :: work_shift, error_psi, error_alph
!
  work_shift = 0.0
! --- Along x-axis
  open(12,file='plot_averaged_error.dat',status='unknown')
  do irg = nrg, 0, -1
    count = 0
    error_psi  = 0.0d0
    error_alph = 0.0d0
    do itg = 0, ntg
      do ipg = 0, npg
        if (irg.gt.irg_exin(itg,ipg).and. &
        &   irg.lt.irg_exout(itg,ipg)) cycle
        error_psi  = error_psi + &
        &     dabs((psi(irg,itg,ipg) - bvxd(irg,itg,ipg))* &
        &     100.0d0/bvxd(irg,itg,ipg))
        error_alph = error_alph + &
        &     dabs((alph(irg,itg,ipg) - bvyd(irg,itg,ipg))* &
        &     100.0d0/bvyd(irg,itg,ipg))
        count = count + 1
      end do
    end do
    write(12,'(1p,10e20.12)') -rg(irg)+work_shift, &
    &     dabs(error_psi)/dble(count), dabs(error_alph)/dble(count)
  end do
  write(12,'(/)')
  do irg = 0, nrg
    count = 0
    error_psi  = 0.0d0
    error_alph = 0.0d0
    do itg = 0, ntg
      do ipg = 0, npg
        if (irg.gt.irg_exin(itg,ipg).and. &
        &   irg.lt.irg_exout(itg,ipg)) cycle
        error_psi  = error_psi + &
        &     dabs((psi(irg,itg,ipg) - bvxd(irg,itg,ipg))* &
        &     100.0d0/bvxd(irg,itg,ipg))
        error_alph = error_alph + &
        &     dabs((alph(irg,itg,ipg) - bvyd(irg,itg,ipg))* &
        &     100.0d0/bvyd(irg,itg,ipg))
        count = count + 1
      end do
    end do
    write(12,'(1p,10e20.12)') rg(irg)+work_shift, &
    &     dabs(error_psi)/dble(count), dabs(error_alph)/dble(count)
  end do
  close(12)
!
end subroutine IO_output_plot_averaged_error
