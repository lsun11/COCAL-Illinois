subroutine sourceterm_surface_int(fnc,irg_surf,sou_surf,dsou_surf)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdr_gridpoint_type0_nosym
!
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_surf(:,:), dsou_surf(:,:)
  real(long), pointer :: fnc_surf(:,:), dfnc_surf(:,:)
  real(long) :: deriv, val
  integer, intent(in) :: irg_surf
  integer :: itg, ipg
!
  call alloc_array2d(fnc_surf, 0, ntg, 0, npg)
  call alloc_array2d(dfnc_surf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      fnc_surf(itg,ipg) = fnc(irg_surf,itg,ipg)
      call grdr_gridpoint_type0_nosym(fnc,deriv,irg_surf,itg,ipg)
      dfnc_surf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,fnc_surf,itg,ipg)
      sou_surf(itg,ipg) = val
      call interpo_linear_type0_2Dsurf(val,dfnc_surf,itg,ipg)
      dsou_surf(itg,ipg) = val
    end do
  end do
  deallocate(fnc_surf)
  deallocate(dfnc_surf)
!
end subroutine sourceterm_surface_int
