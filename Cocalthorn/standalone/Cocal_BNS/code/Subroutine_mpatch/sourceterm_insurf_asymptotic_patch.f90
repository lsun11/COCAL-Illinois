subroutine sourceterm_insurf_asymptotic_patch(impt_bin,impt_apt, &
  &                                           fnc,sou_in,dsou_in)
  use phys_constant, only : long, nnrg
  use grid_parameter, only : nrg, ntg, npg, npgxzm
  use coordinate_grav_extended, only : rgex
!
  use trigonometry_grav_theta, only  : hsinthg
  use trigonometry_grav_phi, only : hcosmpg
!
  use make_array_2d
  use interface_interpo_binary_to_asymptotic_patch
  use interface_interpo_linear_type0_2Dsurf
!
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_in(:,:), dsou_in(:,:)
  real(long), pointer :: fnc_insurf(:,:), dfnc_insurf(:,:)
  integer, intent(IN) :: impt_bin, impt_apt
  real(long) :: deriv, val, r5(5), fr5(5), rv
  real(long) :: rgex_apt(-2:nnrg+2)
  integer :: irg, itg, ipg, ir0, irg0, ii
  integer :: ntg_apt, npg_apt, npgxzm_apt
  real(long), external :: dfdx_4th
!
  call copy_grid_parameter_from_mpt(impt_apt)
  call copy_coordinate_grav_extended_from_mpt(impt_apt)
  ntg_apt = ntg
  npg_apt = npg
  npgxzm_apt = npgxzm
  rgex_apt(-2:nrg+2) = rgex(-2:nrg+2)
  call copy_grid_parameter_from_mpt(impt_bin)
  call copy_coordinate_grav_extended_from_mpt(impt_bin)
!
  call alloc_array2d(fnc_insurf, 0, ntg_apt, 0, npg_apt)
  call alloc_array2d(dfnc_insurf, 0, ntg_apt, 0, npg_apt)
!
  do ipg = 0, npg_apt
    do itg = 0, ntg_apt
      irg = 0
      ir0 = irg - 2
      do ii = 1, 5
        irg0 = ir0 + ii - 1
        call interpo_binary_to_asymptotic_patch(fnc,val,irg0,itg,ipg)
        r5(ii) = rgex_apt(irg0)
        fr5(ii) = val
      end do
      rv = rgex_apt(irg)
      deriv = dfdx_4th(r5,fr5,rv)
      fnc_insurf(itg,ipg) = fr5(3)
      dfnc_insurf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg_apt
    do itg = 1, ntg_apt
! -- to rotate phi by pi, use this
!!      ipgin = ipg + npgxzm_apt  - ((ipg+npgxzm_apt)/(npg_apt+1))*npg_apt
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
end subroutine sourceterm_insurf_asymptotic_patch
