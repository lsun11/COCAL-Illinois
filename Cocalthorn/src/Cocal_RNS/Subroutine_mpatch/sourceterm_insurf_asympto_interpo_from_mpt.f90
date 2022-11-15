subroutine sourceterm_insurf_asympto_interpo_from_mpt &
  &                    (impt_bin,impt_apt,fnc,sou_in,dsou_in)
  use phys_constant, only : long, nnrg
  use grid_parameter, only : nrg, ntg, npg, npgxzm
  use coordinate_grav_extended
  use grid_points_asymptotic_patch, only : ra, tha, phia
  use trigonometry_grav_theta, only  : hsinthg
  use trigonometry_grav_phi, only : hcosmpg
  use make_array_2d
  use interface_interpo_patch_to_active_patch_mpt
  use interface_interpo_linear_type0_2Dsurf
!
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_in(:,:), dsou_in(:,:)
  real(long), pointer :: fnc_insurf(:,:), dfnc_insurf(:,:)
  real(long) :: rc, thc, phic
  integer, intent(in) :: impt_bin, impt_apt
  real(long) :: deriv, val, r5(5), fr5(5), rv
  integer :: irg, itg, ipg, ir0, irg0, ii
  real(long), external :: dfdx_4th
!
  call alloc_array2d(fnc_insurf, 0, ntg, 0, npg)
  call alloc_array2d(dfnc_insurf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      irg = 0
      ir0 = irg - 2
      do ii = 1, 5
        irg0 = ir0 + ii - 1
        rc = ra(irg0,itg,ipg)
        thc = tha(irg0,itg,ipg)
        phic = phia(irg0,itg,ipg)
        call interpo_patch_to_active_patch_mpt(fnc,val,rc,thc,phic)
        r5(ii) = rgex(irg0)
        fr5(ii) = val
      end do
      rv = rgex(irg)
      deriv = dfdx_4th(r5,fr5,rv)
      fnc_insurf(itg,ipg) = fr5(3)
      dfnc_insurf(itg,ipg) = deriv
    end do
  end do
!
! -- ra, tha, phi are rotated by pi for impt = 2 patch
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,fnc_insurf,itg,ipg)
      sou_in(itg,ipg) = val
! sample     sou_in(itg,ipg) = -1.0d0+0.2*hcosmpg(2,ipg)*hsinthg(itg)
      call interpo_linear_type0_2Dsurf(val,dfnc_insurf,itg,ipg)
      dsou_in(itg,ipg) = val
    end do
  end do
!
  deallocate(fnc_insurf)
  deallocate(dfnc_insurf)
end subroutine sourceterm_insurf_asympto_interpo_from_mpt
