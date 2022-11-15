subroutine surf_source_komar_mass(surf_sou_kom,irg)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   ntg, npg
  use def_metric, only  :   tfkijkij, psi, alph
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdr_gridpoint_type0_nosym
  implicit none
  real(long), pointer :: dalph_surf(:,:), surf_sou_kom(:,:)
  integer     ::   irg,itg,ipg
  real(long)  ::   deriv, val
!
  call alloc_array2d(dalph_surf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      call grdr_gridpoint_type0_nosym(alph,deriv,irg,itg,ipg)
      dalph_surf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,dalph_surf,itg,ipg)
      surf_sou_kom(itg,ipg) = val
    end do
  end do
!
  deallocate(dalph_surf)
end subroutine surf_source_komar_mass



