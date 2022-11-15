module def_metric_hij
  implicit none
  real(8), pointer :: hxxd(:,:,:), hxyd(:,:,:), hxzd(:,:,:), &
     &                hyyd(:,:,:), hyzd(:,:,:), hzzd(:,:,:), &
     &                hxxu(:,:,:), hxyu(:,:,:), hxzu(:,:,:), &
     &                hyyu(:,:,:), hyzu(:,:,:), hzzu(:,:,:)
  real(8), pointer :: gaugex(:,:,:), gaugey(:,:,:), gaugez(:,:,:)
end module def_metric_hij
