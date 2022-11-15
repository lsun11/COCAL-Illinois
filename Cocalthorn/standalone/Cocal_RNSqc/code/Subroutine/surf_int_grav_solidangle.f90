subroutine surf_int_grav_solidangle(sousf,surf)
  use phys_constant, only  : long
  use grid_parameter, only : ntg, npg
  use weight_midpoint_grav, only : hwtpgsf
  implicit none
  real(long), pointer     ::  sousf(:,:)
  real(long), intent(out) ::  surf
  integer                 ::  itg, ipg
!
  surf = 0.0d0
  do ipg = 1, npg
    do itg = 1, ntg
      surf = surf + sousf(itg,ipg)*hwtpgsf(itg,ipg)
    end do
  end do
!
end subroutine surf_int_grav_solidangle
