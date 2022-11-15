module def_metric_hij_cartesian
  use phys_constant, only : long
  implicit none
  real(8), pointer :: hxxdca(:,:,:), hxydca(:,:,:), hxzdca(:,:,:), &
     &                hyydca(:,:,:), hyzdca(:,:,:), hzzdca(:,:,:), &
     &                hxxuca(:,:,:), hxyuca(:,:,:), hxzuca(:,:,:), &
     &                hyyuca(:,:,:), hyzuca(:,:,:), hzzuca(:,:,:)
end module def_metric_hij_cartesian
