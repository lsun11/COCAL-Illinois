subroutine rotation_law
  use phys_constant, only  :   long
  use grid_parameter
  use coordinate_grav_r, only : rg
  use def_matter, only : rs, omef, jomef, jomef_int, &
  &                     utf, omeg, jomeg, jomeg_int
  use def_matter_parameter, only : ROT_LAW, ome
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use def_metric_on_SFC_CF, only : alphf, psif, bvxdf, bvydf, bvzdf
  use def_vector_x, only : vec_xf, vec_xg
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long) :: alp_tmp, psi_tmp, by_tmp, hyy_tmp
  real(long) :: omefc, jomefc, jomef_intfc, omegc, jomegc, jomeg_intgc
  real(long) :: xx, yy, zz, Rcyl
  integer    :: irf, itf, ipf
!
  if (ROT_LAW.eq.'DR'.and.nrf.ne.nrf_deform.or. &
  &   ROT_LAW.eq.'OJ'.and.nrf.ne.nrf_deform) then 
    do itf = 0, ntf
      do irf = 0, nrf
        ipf = 0
        xx = vec_xf(irf,itf,ipf,1)
        yy = vec_xf(irf,itf,ipf,2)
        zz = vec_xf(irf,itf,ipf,3)
        Rcyl = sqrt(xx**2 + yy**2)
        omefc = omef(irf,itf,ipf)
!
        alp_tmp = alphf(irf,itf,ipf)
        psi_tmp = psif(irf,itf,ipf)
        by_tmp  = bvydf(irf,itf,ipf)
        hyy_tmp = 0.0d0 
!
        call calc_omega_drot(Rcyl,alp_tmp,psi_tmp,by_tmp,hyy_tmp, &
        &                    omefc,jomefc,jomef_intfc)
        omef(irf,itf,0:npf) = omefc
        jomef(irf,itf,0:npf) = jomefc
        jomef_int(irf,itf,0:npf) = jomef_intfc
!
        xx = vec_xg(irf,itf,ipf,1)
        yy = vec_xg(irf,itf,ipf,2)
        zz = vec_xg(irf,itf,ipf,3)
        Rcyl = sqrt(xx**2 + yy**2)
        omegc = omeg(irf,itf,ipf)
        alp_tmp = alph(irf,itf,ipf)
        psi_tmp = psi(irf,itf,ipf)
        by_tmp  = bvyd(irf,itf,ipf)
        hyy_tmp = 0.0d0
        if (rg(irf).gt.rs(itf,ipf)) then 
           omeg(irf,itf,ipf) = 0.0d0
          jomeg(irf,itf,ipf) = 0.0d0
          jomeg_int(irf,itf,ipf) = 0.0d0
          cycle
        end if
        call calc_omega_drot(Rcyl,alp_tmp,psi_tmp,by_tmp,hyy_tmp, &
        &                    omegc,jomegc,jomeg_intgc)
        omeg(irf,itf,0:npf) = omegc
        jomeg(irf,itf,0:npf) = jomegc
        jomeg_int(irf,itf,0:npf) = jomeg_intgc
      end do
    end do
  else 
    do ipf = 0, npf
      do itf = 0, ntf
        do irf = 0, nrf
          omef(irf,itf,ipf) = ome
         jomef(irf,itf,ipf) = 0.0d0
         jomef_int(irf,itf,ipf) = 0.0d0
          omeg(irf,itf,ipf) = ome
         jomeg(irf,itf,ipf) = 0.0d0
         jomeg_int(irf,itf,ipf) = 0.0d0
         if (rg(irf).gt.rs(itf,ipf)) omeg(irf,itf,ipf) = 0.0d0
        end do
      end do
    end do
  end if
!
end subroutine rotation_law
