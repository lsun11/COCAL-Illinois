subroutine admmass_frg
!
  use phys_constant, only : pi
  use def_metric_1D, only : alphf, psif
  use CB_fR_param_flphy, only : ahoadm
  use def_matter_1D, only : emd, ber, pinx, radi
  use grid_parameter_1D, only : nrf
  use weight_grav_1D, only : wgdr
  implicit none
  real(8) :: alpfc, emdfc, hhfc, prefc, psifc, rhofc, rhoHc, sou, utfc, weiflu
  integer :: ir
!
  ahoadm = 0.0d0
  emd(nrf) = 1.0d-14
  do ir = 0, nrf
!
    if (emd(ir) <= 0.0d0) emd(ir) = 1.0d-14
    emdfc = emd(ir)
    rhofc = emd(ir)**pinx
    prefc = rhofc*emdfc
    alpfc = alphf(ir)
    psifc = psif(ir)
    hhfc  = 1.0d0 + (pinx+1.0d0)*emdfc
    utfc  = hhfc/ber
!
    rhoHc = hhfc*rhofc*(alpfc*utfc)**2 - prefc
    sou = rhoHc*psifc**5
!
    weiflu = wgdr(ir)
    ahoadm = ahoadm + sou*weiflu
!
  end do
!
  ahoadm = radi**3*ahoadm*4.0d0*pi
!
end subroutine admmass_frg
