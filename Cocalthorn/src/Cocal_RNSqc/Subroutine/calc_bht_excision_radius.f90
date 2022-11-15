subroutine calc_bht_excision_radius
  use phys_constant, only : long, pi
  use grid_parameter, only : rgin
  use def_bh_parameter, only : mass_bh, spin_bh
  implicit none
  real(long) :: drh, rh, aom

  rh = mass_bh + dsqrt(mass_bh*mass_bh-spin_bh*spin_bh)
  drh = rh - spin_bh
!  rgin = spin_bh + 0.9*drh

  aom = spin_bh/mass_bh
!  rgin = mass_bh*( 1.0d0 + 0.2d0*dsqrt(1.0d0-aom**2) )
  rgin = mass_bh*( 1.0d0 + 0.8d0*dsqrt(1.0d0-aom**2) )

!  open(120,file="rnspar_rgin.dat", status='unknown')
!  write(120,'(a16,1p,3e23.15)') "r horizon      =", rh 
!  write(120,'(a16,1p,3e23.15)') "Kerr parameter =", spin_bh
!  write(120,'(a16,1p,3e23.15)') "rgin           =", rgin 
!  close(120)
!
end subroutine calc_bht_excision_radius
