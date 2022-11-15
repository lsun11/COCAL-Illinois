subroutine calc_enthalpy_xyzaxis
  use grid_parameter, only : nrf, ntf, npf, ntfeq, ntfpolp, npfxzp, npfyzp
  use def_matter, only : emd
  use def_quantities, only : dhdr_x, dhdr_y, dhdr_z, chi_cusp
  use make_array_3d
  use interface_flgrad_2nd_gridpoint
  implicit none
  real(8), pointer :: hh(:,:,:)
  real(8) :: qq, hhtmp, pre_dum, rho_dum, epsi_dum
  real(8) :: dfdx, dfdy, dfdz, dlnh_x, dlnh_z
  integer :: irf, itf, ipf
!
  call alloc_array3d(hh, 0, nrf, 0, ntf, 0, npf)
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        qq = emd(irf,itf,ipf)
        call peos_q2hprho(qq, hhtmp, pre_dum, rho_dum, epsi_dum)
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
end subroutine calc_enthalpy_xyzaxis
