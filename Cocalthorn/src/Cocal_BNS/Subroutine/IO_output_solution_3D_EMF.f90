subroutine IO_output_solution_3D_EMF
  use phys_constant, only : long
  use def_emfield, only  : va, vaxd, vayd, vazd
  use def_faraday_tensor, only  : fxd_grid, fyd_grid, fzd_grid, fijd_grid
  use grid_parameter, only  :   nrg, ntg, npg
  implicit none
  integer :: irg, itg, ipg
!
! --- EM 1-form and Faraday tensor
  open(13,file='rnsEMF_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        write(13,'(1p,6e20.12)')   va(irg,itg,ipg), &
        &                        vaxd(irg,itg,ipg), &
        &                        vayd(irg,itg,ipg), &
        &                        vazd(irg,itg,ipg)
      end do
    end do
  end do
  close(13)
!
  open(13,file='rnsEMF_faraday_3D.las',status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        write(13,'(1p,6e20.12)') fxd_grid(irg,itg,ipg), &
        &                        fyd_grid(irg,itg,ipg), &
        &                        fzd_grid(irg,itg,ipg), &
        &                        fijd_grid(irg,itg,ipg,1), &
        &                        fijd_grid(irg,itg,ipg,2), &
        &                        fijd_grid(irg,itg,ipg,3)
      end do
    end do
  end do
  close(13)
!
end subroutine IO_output_solution_3D_EMF
