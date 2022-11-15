!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_grav_mpt
  use phys_constant,           only : nnrg, nntg, nnpg, long, nnmpt
  use grid_parameter,          only : nrg, ntg, npg 
  use coordinate_grav_r,       only : drg, rg, hrg
  use coordinate_grav_theta,   only : dthg
  use coordinate_grav_phi,     only : dphig
  use trigonometry_grav_theta, only : sinthg, hsinthg
  use make_array_2d
  use make_array_3d
  implicit none
! weight
  real(long) ::  hwdrg_(nnrg,nnmpt), &
  &              hwdtg_(nntg,nnmpt), &
  &              hwdpg_(nnpg,nnmpt)
  real(long) ::  tzwdrg_(0:nnrg,nnmpt), &
  &              tzwdtg_(0:nntg,nnmpt), &
  &              tzwdpg_(0:nnpg,nnmpt)
  real(long) ::  wdxg_(0:nnrg,nnmpt)
  real(long), pointer :: hwtpgsf_(:,:,:), tzwrtpg_(:,:,:,:), hwrtpg_(:,:,:,:)
end module weight_midpoint_grav_mpt
