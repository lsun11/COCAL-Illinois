subroutine average_error
  use phys_constant, only : long
  use def_metric, only : psi, alph, bvxd, bvyd
  use coordinate_grav_r, only : rg
  use grid_parameter, only  : nrg, ntg, npg
  use grid_points_binary_excision, only : rb
  use grid_parameter_binary_excision, only :ex_radius 
  implicit none
  integer :: irg, itg, ipg, np
  real(long) :: error_psi = 0.0d0, error_alph = 0.0d0
!
  open(12,file='average_error_x.dat',status='unknown')
  do irg = 0, nrg
    error_psi = 0.0d0
    error_alph = 0.0d0
    np = 0
    do itg = 0, ntg
      do ipg = 0, npg
        if (rb(irg,itg,ipg).ge.ex_radius) then
          np = np + 1
          error_psi  = error_psi  + 100*dabs(psi(irg,itg,ipg)-bvxd(irg,itg,ipg))/dabs(bvxd(irg,itg,ipg))     
          error_alph = error_alph + 100*dabs(alph(irg,itg,ipg)-bvyd(irg,itg,ipg))/dabs(bvyd(irg,itg,ipg)) 
        end if    
      end do
    end do
    error_psi = error_psi/dble(np)
    error_alph = error_alph/dble(np)
    write(12,'(1p,6e20.12)')  rg(irg), error_psi, error_alph 
  end do
  close(12)
!
end subroutine average_error
