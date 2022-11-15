subroutine rotation_law_WL
  use phys_constant, only  :   long
  use grid_parameter
  use def_matter, only : rs, omef, jomef, jomef_int, &
  &                     utf, omeg, jomeg, jomeg_int
  use def_matter_parameter, only : ROT_LAW, ome
  use def_metric, only     : alph, psi, bvyd
  use def_metric_hij, only : hyyd
  use def_metric_on_SFC_CF, only : alphf, psif, bvydf
  use def_metric_on_SFC_WL, only : hyydf
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
        hyy_tmp = hyydf(irf,itf,ipf)
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
        hyy_tmp = hyyd(irf,itf,ipf)
        call calc_omega_drot(Rcyl,alp_tmp,psi_tmp,by_tmp,hyy_tmp, &
        &                    omegc,jomegc,jomeg_intgc)
        omeg(irf,itf,0:npf) = omegc
        jomeg(irf,itf,0:npf) = jomegc
        jomeg_int(irf,itf,0:npf) = jomeg_intgc
      end do
    end do
  else 
    omef(0:nrf,0:ntf,0:npf) = ome
    jomef(0:nrf,0:ntf,0:npf) = 0.0d0
    jomef_int(0:nrf,0:ntf,0:npf) = 0.0d0
    omeg(0:nrf,0:ntf,0:npf) = ome
    jomeg(0:nrf,0:ntf,0:npf) = 0.0d0
    jomeg_int(0:nrf,0:ntf,0:npf) = 0.0d0
  end if
!
end subroutine rotation_law_WL
