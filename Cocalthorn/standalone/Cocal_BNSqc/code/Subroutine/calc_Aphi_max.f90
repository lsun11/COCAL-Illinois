subroutine calc_Aphi_max
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf
  use def_emfield, only : vayd
  use def_vector_phi, only : vec_phif
  use integrability_fnc_MHD, only : Aphi_max_surf, Aphi_max_vol
  use interface_interpo_gr2fl_type0
  implicit none
  real(long) :: Aphi_tmp, Ay
  integer    :: irf, itf, ipf
!
  Aphi_max_surf = 0.0d0
  irf = nrf; ipf = 0
  do itf = 0, ntf
    call interpo_gr2fl_type0(Ay,vayd,irf,itf,ipf)
    Aphi_tmp      = Ay*vec_phif(irf,itf,ipf,2)
    Aphi_max_surf = dmax1(Aphi_tmp,Aphi_max_surf)
  end do
!
  Aphi_max_vol = 0.0d0
  ipf = 0
  do itf = 0, ntf
    do irf = 0, nrf
      call interpo_gr2fl_type0(Ay,vayd,irf,itf,ipf)
      Aphi_tmp = Ay*vec_phif(irf,itf,ipf,2)
      Aphi_max_vol = dmax1(Aphi_tmp,Aphi_max_vol)
    end do
  end do
!
end subroutine calc_Aphi_max
