subroutine test_excurve
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use def_metric, only : psi, tfkij
  use def_metric_excurve_grid, only : tfkij_grid
  use make_array_2d
  implicit none
  real(long), pointer :: psi_bhsurf(:,:), dpsi_bhsurf(:,:)
  real(long) :: work(2,2), val
  integer    :: irg, itg, ipg, ii, jj, kk, irg1, irg2, ipg1, ipg2, itg1, itg2
!
  open(15,file='kij_r_line.txt',status='unknown')
  ipg1 = npg       !  refers to mid points for tfkij
  ipg2 = 1         !  refers to mid points for tfkij
  itg1 = ntg/2     !  refers to mid points for tfkij
  itg2 = ntg/2+1   !  refers to mid points for tfkij
  ipg  = 0         !  refers to grid points for tfkij_grid
  itg  =ntg/2      !  refers to grid points for tfkij_grid

  do ii = 1,3
    do jj = ii,3
      write(15,'(a4,i5,a10,i5)')  '# i=',ii, '        j=',jj      
      write(15,'(1p,2e20.12)')  rg(0), tfkij_grid(0,itg,ipg,ii,jj)
      do irg = 1,nrg-1
        val = 0.125d0*(tfkij(irg,itg1,ipg1,ii,jj) + tfkij(irg,itg2,ipg1,ii,jj) + &
             &         tfkij(irg,itg1,ipg2,ii,jj) + tfkij(irg,itg2,ipg2,ii,jj) + &
             &         tfkij(irg+1,itg1,ipg1,ii,jj) + tfkij(irg+1,itg2,ipg1,ii,jj) + &
             &         tfkij(irg+1,itg1,ipg2,ii,jj) + tfkij(irg+1,itg2,ipg2,ii,jj) )

        write(15,'(1p,3e20.12)')  rg(irg), tfkij_grid(irg,itg,ipg,ii,jj), val
      end do
      write(15,'(1p,2e20.12)')  rg(nrg), tfkij_grid(nrg,itg,ipg,ii,jj)
      write(15,'(a1)')  ' '
    end do
  end do
  close(15)

!
!
end subroutine test_excurve
