subroutine calc_enthalpy_xyzaxis_qeos
  use grid_parameter, only : nrf, ntf, npf, ntfeq, ntfpolp, npfxzp, npfyzp
  use def_matter, only : rhof
  use def_quantities, only : dhdr_x, dhdr_y, dhdr_z, chi_cusp
  use make_array_3d
  use interface_flgrad_2nd_gridpoint
  implicit none
  real(8), pointer :: hh(:,:,:)
  real(8) :: hhtmp, pre_c, rho_c, epsi_c, dummy
  real(8) :: dfdx, dfdy, dfdz, dlnh_x, dlnh_z
  integer :: irf, itf, ipf
!
  call alloc_array3d(hh, 0, nrf, 0, ntf, 0, npf)
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        rho_c = rhof(irf,itf,ipf)
        call quark_rho2phenedpdrho(rho_c, pre_c, hhtmp, epsi_c, dummy)
        hh(irf,itf,ipf) = hhtmp
      end do
    end do
  end do
!
  call flgrad_2nd_gridpoint(hh,dfdx,dfdy,dfdz,nrf,ntfeq,npfxzp)
  dhdr_x = dfdx
  call flgrad_2nd_gridpoint(hh,dfdx,dfdy,dfdz,nrf,ntfeq,npfyzp)
  dhdr_y = dfdy
  call flgrad_2nd_gridpoint(hh,dfdx,dfdy,dfdz,nrf,ntfpolp,0)
  dhdr_z = dfdz
!
  dlnh_x   = dhdr_x/hh(nrf,ntfeq,npfxzp)
  dlnh_z   = dhdr_z/hh(nrf,ntfpolp,0)
  chi_cusp = dlnh_x/dlnh_z
!
  deallocate(hh)
!
end subroutine calc_enthalpy_xyzaxis_qeos
