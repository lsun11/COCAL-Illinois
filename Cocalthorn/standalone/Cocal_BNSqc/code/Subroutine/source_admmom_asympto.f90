subroutine source_admmom_asympto(sousfv,irg)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  use def_metric_excurve_grid, only : tfkij_grid
  use coordinate_grav_r, only : rg
  use def_vector_x, only   : vec_xg
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sousfv(:,:,:)
  integer    :: irg
  real(long), pointer :: fnc(:,:)
  real(long) :: na, eb, Aij, val, eia(3,3)
  integer    :: itg, ipg, ia, ib, ii
!
  call alloc_array2d(fnc, 0,ntg, 0,npg)
  eia(1,1) = 1.0d0 ; eia(1,2) = 0.0d0 ; eia(1,3) = 0.0d0
  eia(2,1) = 0.0d0 ; eia(2,2) = 1.0d0 ; eia(2,3) = 0.0d0
  eia(3,1) = 0.0d0 ; eia(3,2) = 0.0d0 ; eia(3,3) = 1.0d0
!
do ii = 1, 3
!
  do ipg = 0, npg
    do itg = 0, ntg
      fnc(itg,ipg) = 0.0d0
      do ib = 1, 3
        do ia = 1, 3
          Aij = tfkij_grid(irg,itg,ipg,ia,ib)
          eb = eia(ii,ib)
          na = vec_xg(irg,itg,ipg,ia)/rg(irg)
          fnc(itg,ipg) = fnc(itg,ipg) + Aij*eb*na
        end do
      end do
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,fnc,itg,ipg)
      sousfv(itg,ipg,ii) = val
    end do
  end do
!
end do
!
  deallocate(fnc)
!
end subroutine source_admmom_asympto
