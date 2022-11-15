module def_Lie_derivatives
  implicit none
  real(8), pointer :: elpxx(:,:,:), elpxy(:,:,:), elpxz(:,:,:), &
  &                   elpyy(:,:,:), elpyz(:,:,:), elpzz(:,:,:)
  real(8), pointer :: rlpxx(:,:,:), rlpxy(:,:,:), rlpxz(:,:,:), &
  &                   rlpyy(:,:,:), rlpyz(:,:,:), rlpzz(:,:,:)
  real(8), pointer :: rlbxx(:,:,:), rlbxy(:,:,:), rlbxz(:,:,:), &
  &                   rlbyy(:,:,:), rlbyz(:,:,:), rlbzz(:,:,:)
end module def_Lie_derivatives
