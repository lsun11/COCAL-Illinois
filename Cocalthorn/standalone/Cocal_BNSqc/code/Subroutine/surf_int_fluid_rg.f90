subroutine surf_int_fluid_rg(souf,surf,irf)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrf, ntf, npf
  use weight_midpoint_fluid, only : hwtpfsf
  use def_matter, only : rs
  use coordinate_grav_r, only : rg
  use interface_interpo_linear_type0
  use interface_interpo_linear_type0_2Dsurf
  implicit none
!
  integer     ::   irf,itf,ipf
  real(long)  ::   surf, hsouf, hsurf
  real(long),pointer  ::   souf(:,:)
!
  surf = 0.0d0
  do ipf = 1, npf
    do itf = 1, ntf
!      call interpo_linear_type0(hsouf,souf,irf,itf,ipf)
      hsouf = souf(itf,ipf)
      call interpo_linear_type0_2Dsurf(hsurf,rs,itf,ipf)
      surf = surf + hsouf*hsurf**2*hwtpfsf(itf,ipf)
    end do
  end do
  surf = surf*rg(irf)**2
!
end subroutine surf_int_fluid_rg
