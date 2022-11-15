subroutine surf_int_grav(sousf,surf,irg)
  use phys_constant, only  : long
  use grid_parameter, only    : ntg, npg
  use coordinate_grav_r, only : hrg
  use weight_midpoint_grav, only : hwtpgsf
  implicit none
  real(long), pointer     ::  sousf(:,:)
  real(long), intent(out) ::  surf
  integer, intent(in)     ::  irg
  integer                 ::  itg, ipg
  surf = 0.0d0
  do ipg = 1, npg
    do itg = 1, ntg
      surf = surf + sousf(itg,ipg)*hwtpgsf(itg,ipg)
    end do
  end do
  surf = surf*hrg(irg)**2
end subroutine surf_int_grav
