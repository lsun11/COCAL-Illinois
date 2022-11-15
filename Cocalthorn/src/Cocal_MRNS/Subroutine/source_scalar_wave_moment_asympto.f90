subroutine source_scalar_wave_moment_asympto(sousfv,irg)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  use def_metric, only : psi
  use coordinate_grav_r, only : rg
  use def_vector_x, only   : vec_xg
  use make_array_2d
  use interface_grgrad_gridpoint_4th_type0
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sousfv(:,:,:)
  integer    :: irg
  real(long), pointer :: fnc(:,:)
  real(long) :: na, psigc, val, dfdx, dfdy, dfdz
  integer    :: itg, ipg, ia, ib, ii
!
  call alloc_array2d(fnc, 0,ntg, 0,npg)
!
do ii = 1, 3
!
  do ipg = 0, npg
    do itg = 0, ntg
!xpsi      psigc = psi(irg,itg,ipg)
!xpsi      na = vec_xg(irg,itg,ipg,ii)/rg(irg)
!xpsi      fnc(itg,ipg) = psigc*na
      call grgrad_gridpoint_4th_type0(psi,dfdx,dfdy,dfdz,irg,itg,ipg,'ns')
      if (ii.eq.1) fnc(itg,ipg) = dfdx
      if (ii.eq.2) fnc(itg,ipg) = dfdy
      if (ii.eq.3) fnc(itg,ipg) = dfdz
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
end subroutine source_scalar_wave_moment_asympto
