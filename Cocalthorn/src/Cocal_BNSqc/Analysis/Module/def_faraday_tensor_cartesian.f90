module def_faraday_tensor_cartesian
  use phys_constant, only : long
  implicit none
  real(long), pointer :: fxdca(:,:,:), fydca(:,:,:), fzdca(:,:,:)
  real(long), pointer :: fxuca(:,:,:), fyuca(:,:,:), fzuca(:,:,:)
  real(long), pointer :: fijdca(:,:,:,:),fijuca(:,:,:,:)
  real(long), pointer :: fijduca(:,:,:,:,:),fidfiuca(:,:,:),fijfijca(:,:,:)
  real(long), pointer :: fxd_gridca(:,:,:),fyd_gridca(:,:,:),fzd_gridca(:,:,:)
  real(long), pointer :: fxu_gridca(:,:,:),fyu_gridca(:,:,:),fzu_gridca(:,:,:)
  real(long), pointer :: fijd_gridca(:,:,:,:), fiju_gridca(:,:,:,:)
  real(long), pointer :: Lie_bFxdca(:,:,:),Lie_bFydca(:,:,:),Lie_bFzdca(:,:,:)
end module def_faraday_tensor_cartesian
