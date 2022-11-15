module def_matter_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  emdg_(:,:,:,:), utg_(:,:,:,:)
  real(long), pointer  ::  emd_(:,:,:,:),  utf_(:,:,:,:)
  real(long), pointer  ::  rs_(:,:,:),  lambda_(:,:,:,:)
  real(long), pointer  ::  rhof_(:,:,:,:)
  real(long), pointer  ::  uxf_(:,:,:,:), uyf_(:,:,:,:), uzf_(:,:,:,:)
  real(long), pointer  ::  uxg_(:,:,:,:), uyg_(:,:,:,:), uzg_(:,:,:,:)
  
  real(long), pointer  ::  omeg_(:,:,:,:), jomeg_(:,:,:,:), jomeg_int_(:,:,:,:)
  real(long), pointer  ::  omef_(:,:,:,:), jomef_(:,:,:,:), jomef_int_(:,:,:,:)

  real(long), pointer  ::  vep_(:,:,:,:)
  real(long), pointer  ::  vepxf_(:,:,:,:), vepyf_(:,:,:,:), vepzf_(:,:,:,:)
  real(long), pointer  ::  vepxg_(:,:,:,:), vepyg_(:,:,:,:), vepzg_(:,:,:,:)
end module def_matter_mpt
