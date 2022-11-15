subroutine sourceterm_outsurf_interpo_from_asympto_mpt &
  &            (impt_bin,impt_apt,fnc,sou_out,dsou_out)
  use phys_constant, only : long, nnrg
  use grid_parameter, only : nrg, ntg, npg, npgxzm
  use coordinate_grav_extended
  use grid_points_binary_in_asympto, only : rb_a, thb_a, phib_a
  use trigonometry_grav_theta, only  : hsinthg
  use trigonometry_grav_phi, only : hcosmpg
  use make_array_2d
  use interface_interpo_patch_to_active_patch_mpt
  use interface_interpo_linear_type0_2Dsurf
!
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_out(:,:), dsou_out(:,:)
  real(long), pointer :: fnc_outsurf(:,:), dfnc_outsurf(:,:)
  integer, intent(in) :: impt_bin, impt_apt
  real(long) :: deriv, val, r5(5), fr5(5), rv, rc, thc, phic
  integer :: irg, itg, ipg, ir0, irg0, ii
  real(long), external :: dfdx_4th
!
  call alloc_array2d(fnc_outsurf, 0, ntg, 0, npg)
  call alloc_array2d(dfnc_outsurf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      irg = nrg
      ir0 = irg - 2
      do ii = 1, 5
        irg0 = ir0 + ii - 1
        rc = rb_a(irg0,itg,ipg)
        thc = thb_a(irg0,itg,ipg)
        phic = phib_a(irg0,itg,ipg)
        call interpo_patch_to_active_patch_mpt(fnc,val,rc,thc,phic)
        r5(ii) = rgex(irg0)
        fr5(ii) = val
      end do
      rv = rgex(irg)
      deriv = dfdx_4th(r5,fr5,rv)
      fnc_outsurf(itg,ipg) = fr5(3)
      dfnc_outsurf(itg,ipg) = deriv
    end do
  end do
!
! -- phib_a is rotated by pi for impt = 2 patch
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,fnc_outsurf,itg,ipg)
      sou_out(itg,ipg) = val
! sample     sou_out(itg,ipg) = -1.0d0+0.2*hcosmpg(2,ipg)*hsinthg(itg)
      call interpo_linear_type0_2Dsurf(val,dfnc_outsurf,itg,ipg)
      dsou_out(itg,ipg) = val
    end do
  end do
!
  deallocate(fnc_outsurf)
  deallocate(dfnc_outsurf)
end subroutine sourceterm_outsurf_interpo_from_asympto_mpt
