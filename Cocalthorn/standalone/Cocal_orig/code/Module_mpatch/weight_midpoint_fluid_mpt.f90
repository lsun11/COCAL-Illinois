!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_fluid_mpt
  use phys_constant,           only : nnrf, nntf, nnpf, long, nnmpt
  use grid_parameter,          only : nrf, ntf, npf, nrg, ntg, npg
  use coordinate_grav_r,       only : drg, rg, hrg
  use coordinate_grav_theta,   only : dthg
  use coordinate_grav_phi,     only : dphig
  use trigonometry_grav_theta, only : hsinthg, sinthg
  use make_array_4d
  implicit none
! weight
  Real(long)          ::  hwdrf_(nnrf,nnmpt), hwdtf_(nntf,nnmpt), hwdpf_(nnpf,nnmpt)
  Real(long)          ::  tzwdrf_(0:nnrf,nnmpt), siwdrf_(0:nnrf,nnmpt)
  Real(long)          ::  siwdtf_(0:nntf,nnmpt), tzwdpf_(0:nnpf,nnmpt)
  Real(long)          ::  wdxf_(0:nnrf,nnmpt)
  Real(long), pointer ::  hwrtpf_(:,:,:,:), tzwrtpf_(:,:,:,:)
  Real(long), pointer ::  siwrtpf_(:,:,:,:), rtsiwrtpf_(:,:,:,:)
end module weight_midpoint_fluid_mpt
