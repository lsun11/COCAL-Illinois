subroutine source_angmom_asympto(sousf,irg)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  use def_metric_excurve_grid, only : tfkij_grid
  use coordinate_grav_r, only : rg
  use def_vector_x, only   : vec_xg
  use def_vector_phi, only : vec_phig
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sousf(:,:)
  integer    :: irg
  real(long), pointer :: fnc(:,:)
  real(long) :: na, vphi_cm, Aij, val
  integer    :: itg, ipg, ia, ib
!
  call alloc_array2d(fnc, 0,ntg, 0,npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      fnc(itg,ipg) = 0.0d0
      do ib = 1, 3
        do ia = 1, 3
          Aij = tfkij_grid(irg,itg,ipg,ia,ib)
          vphi_cm = vec_phig(irg,itg,ipg,ib)
          na = vec_xg(irg,itg,ipg,ia)/rg(irg)
          fnc(itg,ipg) = fnc(itg,ipg) + Aij*vphi_cm*na
        end do
      end do
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,fnc,itg,ipg)
      sousf(itg,ipg) = val
    end do
  end do
!
  deallocate(fnc)
!
end subroutine source_angmom_asympto
