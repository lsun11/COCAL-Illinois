subroutine source_charge_asympto(sousf,irg)
  use phys_constant, only  : long, pi
  use grid_parameter, only : ntg, npg
  use coordinate_grav_r, only: rg
  use def_metric,     only : psi
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_faraday_tensor, only : fxd_grid, fyd_grid, fzd_grid
  use def_vector_x, only : vec_xg
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: dfnc(:,:), sousf(:,:)
  integer    :: irg, itg, ipg
  real(long) :: deriv, val 
  real(long) :: psi2, fxdc, fydc, fzdc, erx, ery, erz
  real(long) :: gmxxu, gmxyu, gmxzu, gmyxu, gmyyu, gmyzu, gmzxu, gmzyu, gmzzu
!
  call alloc_array2d(dfnc, 0,ntg, 0,npg)
  call calc_vector_x_grav(1)
!
  do ipg = 0, npg
    do itg = 0, ntg
      psi2 = psi(irg,itg,ipg)**2
      fxdc = fxd_grid(irg,itg,ipg)
      fydc = fyd_grid(irg,itg,ipg)
      fzdc = fzd_grid(irg,itg,ipg)
      erx  = vec_xg(irg,itg,ipg,1)/rg(irg)
      ery  = vec_xg(irg,itg,ipg,2)/rg(irg)
      erz  = vec_xg(irg,itg,ipg,3)/rg(irg)
      gmxxu= hxxu(irg,itg,ipg) + 1.0d0
      gmxyu= hxyu(irg,itg,ipg)
      gmxzu= hxzu(irg,itg,ipg)
      gmyyu= hyyu(irg,itg,ipg) + 1.0d0
      gmyzu= hyzu(irg,itg,ipg)
      gmzzu= hzzu(irg,itg,ipg) + 1.0d0
      gmyxu= gmxyu
      gmzxu= gmxzu
      gmzyu= gmyzu
!
      dfnc(itg,ipg) = (gmxxu*fxdc*erx + gmxyu*fxdc*ery + gmxzu*fxdc*erz &
      &              + gmyxu*fydc*erx + gmyyu*fydc*ery + gmyzu*fydc*erz &
      &              + gmzxu*fzdc*erx + gmzyu*fzdc*ery + gmzzu*fzdc*erz)*psi2
!
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,dfnc,itg,ipg)
      sousf(itg,ipg) = val
    end do
  end do
!
  deallocate(dfnc)
end subroutine source_charge_asympto
