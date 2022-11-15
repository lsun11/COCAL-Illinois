subroutine source_fnc_moment_asympto(fnc,sousfv,irg)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  use coordinate_grav_r, only : rg
  use def_vector_x, only   : vec_xg
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: fnc(:,:,:), sousfv(:,:,:)
  real(long), pointer :: sousf(:,:)
  integer    :: irg
  real(long) :: na, fncgc, val, dfdx, dfdy, dfdz
  integer    :: itg, ipg, ia, ib, ii
!
  call alloc_array2d(sousf, 0,ntg, 0,npg)
!
  do ii = 1, 3
!
    do ipg = 0, npg
      do itg = 0, ntg
        fncgc = fnc(irg,itg,ipg)
!        na = vec_xg(irg,itg,ipg,ii)/rg(irg)
        na = vec_xg(irg,itg,ipg,ii)
        sousf(itg,ipg) = fncgc*na
      end do
    end do
!
    do ipg = 1, npg
      do itg = 1, ntg
        call interpo_linear_type0_2Dsurf(val,sousf,itg,ipg)
        sousfv(itg,ipg,ii) = val
      end do
    end do
!
  end do
!
  deallocate(sousf)
!
end subroutine source_fnc_moment_asympto
