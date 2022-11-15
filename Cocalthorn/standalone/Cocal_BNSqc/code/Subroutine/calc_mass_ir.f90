subroutine calc_mass_ir(mass_ir)
  use phys_constant, only  : long
  use grid_parameter, only    : nrg
  use coordinate_grav_r, only : rg
  implicit none
  integer   ::  irg, mass_ir
  real(long) :: mass_r
!
  mass_r = 1.0d-02*rg(nrg)
!
  do irg = 0, nrg-1
    if(mass_r.ge.rg(irg).and.mass_r.lt.rg(irg+1)) then
      mass_ir = irg
      exit
    end if
  end do
!
end subroutine calc_mass_ir
