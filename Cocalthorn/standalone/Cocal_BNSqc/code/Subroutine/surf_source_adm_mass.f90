subroutine surf_source_adm_mass(surf_sou_adm,irg)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   ntg, npg
  use def_metric, only  :   tfkijkij, psi, alph
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdr_gridpoint_type0_nosym
  implicit none
  real(long), pointer :: dpsi_surf(:,:), surf_sou_adm(:,:)
  integer     ::   irg,itg,ipg
  real(long)  ::   deriv, val 
!
  call alloc_array2d(dpsi_surf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      call grdr_gridpoint_type0_nosym(psi,deriv,irg,itg,ipg)
      dpsi_surf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,dpsi_surf,itg,ipg)
      surf_sou_adm(itg,ipg) = val
    end do
  end do
!
  deallocate(dpsi_surf)
end subroutine surf_source_adm_mass



