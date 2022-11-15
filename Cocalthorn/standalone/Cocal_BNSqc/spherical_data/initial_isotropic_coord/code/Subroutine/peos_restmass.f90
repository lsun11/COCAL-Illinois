subroutine peos_restmass
!
  use phys_constant, only : pi
  use def_metric_1D, only : alphf, psif
  use CB_fR_param_flphy, only : ahores
  use def_matter_1D, only : emd, ber, pinx, radi
  use grid_parameter_1D, only : nrf
  use weight_grav_1D, only : wgdr
  implicit none
  real(8) :: alpfc, emdfc, hhfc, psifc, utfc, weiflu, rhofc, ene, prefc
  integer :: ir
!
  ahores = 0.0d0
  emd(nrf) = 1.0d-14
  do ir = 0, nrf
!
    if (emd(ir) <= 0.0d0) emd(ir) = 1.0d-14
    emdfc = emd(ir)
!    rhofc = emd(ir)**pinx
    alpfc = alphf(ir)
    psifc = psif(ir)
!    hhfc  = 1.0d0 + (pinx+1.0d0)*emdfc
    call peos_q2hprho(emdfc, hhfc, prefc, rhofc, ene)
    utfc  = hhfc/ber
!
    weiflu = wgdr(ir)
    ahores = ahores + rhofc*alpfc*utfc*psifc**6*weiflu
!
  end do
!
  ahores = radi**3*ahores*4.0d0*pi
!
end subroutine peos_restmass
