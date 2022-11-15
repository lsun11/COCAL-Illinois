subroutine source_AHarea_CF(sousf)
  use phys_constant, only  : long, pi
  use grid_parameter, only : ntg, npg
  use def_metric, only  : psi
  use interface_interpo_linear_surface_type0
  implicit none
  real(long), pointer  :: sousf(:,:)
  real(long) :: val
  integer :: irg, itg, ipg
!
  irg = 0
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_surface_type0(val,psi,irg,itg,ipg)
      sousf(itg,ipg) = val**4
    end do
  end do
!
end subroutine source_AHarea_CF
