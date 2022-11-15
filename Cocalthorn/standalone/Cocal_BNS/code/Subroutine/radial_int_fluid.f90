subroutine radial_int_fluid(sou,radius,it,ip)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrf
  use coordinate_grav_r, only : drg
  use weight_midpoint_fluid, only : wdxf
  use def_matter, only : rs
  implicit none
!
  real(long)  :: souf, surf
  integer     :: ir
  real(long), pointer     :: sou(:)
  real(long), intent(out) :: radius
  integer,    intent(in)  :: it, ip
!
! -- Mid-point integration overestimates the mass by 1% for nr = 16.
!
  radius = 0.0d0
  do ir = 0, nrf
    souf = sou(ir)
    surf = rs(it,ip)
    radius = radius + souf * surf * wdxf(ir)
  end do
!
end subroutine radial_int_fluid
