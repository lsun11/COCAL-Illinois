subroutine source_qua_loc_spin_peos(soug,irs)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf
  use make_array_2d
  use def_metric, only :  psi, tfkij
  use coordinate_grav_r, only : hrg
  use def_vector_x, only : hvec_xg
  use def_vector_phi, only : hvec_phig
  use interface_interpo_linear_type0
  use make_array_2d
  implicit none
  real(long), pointer :: soug(:,:)
  real(long) :: ni, vphi_bh, vphi_cm, Aij_surf, psi6
  real(long) :: x(2),y(2), val
  integer    :: irg, itg, ipg, ia, ib, irs, ic
!
!  irs is a midpoint
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0(val,psi,irs,itg,ipg)
      psi6 = val**6
      soug(itg,ipg)=0.0d0
      do ib = 1, 3
        do ia = 1, 3
          Aij_surf = tfkij(irs,itg,ipg,ia,ib)
          vphi_cm = hvec_phig(irs,itg,ipg,ib)
          ni = +hvec_xg(irs,itg,ipg,ia)/hrg(irs)
          soug(itg,ipg) = soug(itg,ipg) + Aij_surf*vphi_cm*ni
        end do
      end do
!      if (irs==(nrf+1).and.itg==ntg.and.ipg==1) then
!        write(6,'(a7,1p,7e14.6)') "sou,psi",  soug(itg,ipg), psi6
!        write(6,'(a7,1p,7e14.6)') "tv,nv= ",  &
!           & hvec_phig(irs,itg,ipg,1), hvec_phig(irs,itg,ipg,2), hvec_phig(irs,itg,ipg,3), &
!           & hvec_xg(irs,itg,ipg,1)/hrg(irs), hvec_xg(irs,itg,ipg,2)/hrg(irs), hvec_xg(irs,itg,ipg,3)/hrg(irs)
!        write(6,'(a7,1p,7e14.6)') "tfkij= ", tfkij(irs,itg,ipg,1,1), tfkij(irs,itg,ipg,1,2), tfkij(irs,itg,ipg,1,3), &
!         &                     tfkij(irs,itg,ipg,2,2), tfkij(irs,itg,ipg,2,3), tfkij(irs,itg,ipg,3,3)
!      end if

      soug(itg,ipg) = soug(itg,ipg)*psi6
    end do
  end do
!
end subroutine source_qua_loc_spin_peos

