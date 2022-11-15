module def_faraday_tensor
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  fxd(:,:,:), fyd(:,:,:), fzd(:,:,:)
  real(long), pointer  ::  fxu(:,:,:), fyu(:,:,:), fzu(:,:,:)
  real(long), pointer  ::  fijd(:,:,:,:), fiju(:,:,:,:)
  real(long), pointer  ::  fijdu(:,:,:,:,:), fidfiu(:,:,:), fijfij(:,:,:)
  real(long), pointer  ::  fxd_grid(:,:,:), fyd_grid(:,:,:), fzd_grid(:,:,:)
  real(long), pointer  ::  fxu_grid(:,:,:), fyu_grid(:,:,:), fzu_grid(:,:,:)
  real(long), pointer  ::  fijd_grid(:,:,:,:), fiju_grid(:,:,:,:)
  real(long), pointer  ::  fijdu_grid(:,:,:,:,:), fidfiu_grid(:,:,:)
  real(long), pointer  ::  fijfij_grid(:,:,:)
  real(long), pointer  ::  Lie_bFxd(:,:,:), Lie_bFyd(:,:,:), Lie_bFzd(:,:,:)
end module def_faraday_tensor
