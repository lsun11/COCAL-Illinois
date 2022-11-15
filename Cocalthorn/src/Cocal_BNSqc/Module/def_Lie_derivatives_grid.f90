module def_Lie_derivatives_grid
  implicit none
  real(8), pointer :: elpxx_grid(:,:,:), elpxy_grid(:,:,:), &
  &                   elpxz_grid(:,:,:), elpyy_grid(:,:,:), &
  &                   elpyz_grid(:,:,:), elpzz_grid(:,:,:)
  real(8), pointer :: rlpxx_grid(:,:,:), rlpxy_grid(:,:,:), &
  &                   rlpxz_grid(:,:,:), rlpyy_grid(:,:,:), &
  &                   rlpyz_grid(:,:,:), rlpzz_grid(:,:,:)
  real(8), pointer :: rlbxx_grid(:,:,:), rlbxy_grid(:,:,:), &
  &                   rlbxz_grid(:,:,:), rlbyy_grid(:,:,:), &
  &                   rlbyz_grid(:,:,:), rlbzz_grid(:,:,:)
end module def_Lie_derivatives_grid
