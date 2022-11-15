! par takes values +1 or -1
subroutine sourceterm_exsurf_eqm_binary_parity(fnc,sou_ex,dsou_ex,par)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg, npgxzm
  use grid_parameter_binary_excision, only : ex_nrg
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdr_gridpoint_type0
!
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_ex(:,:), dsou_ex(:,:)
  real(long), pointer :: fnc_exsurf(:,:), dfnc_exsurf(:,:)
  real(long) :: deriv, val, par
  integer :: itg, ipg, itgex, ipgex
!
  call alloc_array2d(fnc_exsurf, 0, ntg, 0, npg)
  call alloc_array2d(dfnc_exsurf, 0, ntg, 0, npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      fnc_exsurf(itg,ipg) = fnc(ex_nrg,itg,ipg)
      call grdr_gridpoint_type0(fnc,deriv,ex_nrg,itg,ipg)
      dfnc_exsurf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
!
      itgex = itg
      ipgex = ipg + npgxzm  - ((ipg+npgxzm)/(npg+1))*npg
!reflection sym      ipgex = npgxzm + 1 - ipg + ((ipg+npgxzm)/(npg+1))*npg
!
!!test      call interpo_linear_type0_2Dsurf(val,fnc_exsurf,itgex,ipgex)
!!test      sou_ex(itg,ipg) = val
!!test      call interpo_linear_type0_2Dsurf(val,dfnc_exsurf,itgex,ipgex)
!!test      dsou_ex(itg,ipg) = val
      call interpo_linear_type0_2Dsurf(val,fnc_exsurf,itg,ipg)
      sou_ex(itgex,ipgex) = par*val
      call interpo_linear_type0_2Dsurf(val,dfnc_exsurf,itg,ipg)
      dsou_ex(itgex,ipgex) = par*val
    end do
  end do
  deallocate(fnc_exsurf)
  deallocate(dfnc_exsurf)
!
end subroutine sourceterm_exsurf_eqm_binary_parity
