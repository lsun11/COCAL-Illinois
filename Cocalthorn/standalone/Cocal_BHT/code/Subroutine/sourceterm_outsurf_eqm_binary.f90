subroutine sourceterm_outsurf_eqm_binary(fnc,sou_out,dsou_out)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, npgxzm
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdr_gridpoint_type0
!
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_out(:,:), dsou_out(:,:)
  real(long), pointer :: fnc_outsurf(:,:), dfnc_outsurf(:,:)
  real(long) :: deriv, val
  integer :: itg, ipg, itgout, ipgout
!
  call alloc_array2d(fnc_outsurf, 0, ntg, 0, npg)
  call alloc_array2d(dfnc_outsurf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      fnc_outsurf(itg,ipg) = fnc(nrg,itg,ipg)
      call grdr_gridpoint_type0(fnc,deriv,nrg,itg,ipg)
      dfnc_outsurf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
!
      itgout = itg
      ipgout = ipg + npgxzm  - ((ipg+npgxzm)/(npg+1))*npg
      call interpo_linear_type0_2Dsurf(val,fnc_outsurf,itg,ipg)
      sou_out(itgout,ipgout) = val
      call interpo_linear_type0_2Dsurf(val,dfnc_outsurf,itg,ipg)
      dsou_out(itgout,ipgout) = val
    end do
  end do
  deallocate(fnc_outsurf)
  deallocate(dfnc_outsurf)
!
end subroutine sourceterm_outsurf_eqm_binary
