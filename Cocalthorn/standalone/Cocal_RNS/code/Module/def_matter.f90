module def_matter
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  emd(:,:,:),  utf(:,:,:), omef(:,:,:)
  real(long), pointer  ::  emdg(:,:,:), utg(:,:,:), omeg(:,:,:)
  real(long), pointer  ::  jomef(:,:,:), jomef_int(:,:,:)
  real(long), pointer  ::  jomeg(:,:,:), jomeg_int(:,:,:)
  real(long), pointer  ::  rs(:,:),  lambda(:,:,:)
  real(long), pointer  ::  ergoin(:,:), ergoout(:,:)
  real(long), pointer  ::  rhog(:,:,:)
  real(long), pointer  ::  rhof(:,:,:), hhf(:,:,:)
  real(long), pointer  ::  utdf(:,:,:), utdg(:,:,:)
  real(long), pointer  ::  uxf(:,:,:),  uyf(:,:,:),  uzf(:,:,:)
  real(long), pointer  ::  uxdf(:,:,:), uydf(:,:,:), uzdf(:,:,:)
  real(long), pointer  ::  uxg(:,:,:),  uyg(:,:,:),  uzg(:,:,:)
  real(long), pointer  ::  uxdg(:,:,:), uydg(:,:,:), uzdg(:,:,:)
  real(long), pointer  ::  vep(:,:,:)
  real(long), pointer  ::  vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
  real(long), pointer  ::  vepxg(:,:,:), vepyg(:,:,:), vepzg(:,:,:)
end module def_matter

