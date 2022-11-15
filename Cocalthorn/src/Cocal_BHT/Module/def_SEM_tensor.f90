module def_SEM_tensor
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  rhoH(:,:,:)
  real(long), pointer  ::  jmd(:,:,:,:), jmu(:,:,:,:)
  real(long), pointer  ::  smijd(:,:,:,:,:), smiju(:,:,:,:,:)
  real(long), pointer  ::  trsm(:,:,:)
  real(long), pointer  ::  rhoH_grid(:,:,:)
  real(long), pointer  ::  jmd_grid(:,:,:,:), jmu_grid(:,:,:,:)
  real(long), pointer  ::  smijd_grid(:,:,:,:,:), smiju_grid(:,:,:,:,:)
  real(long), pointer  ::  trsm_grid(:,:,:)
end module def_SEM_tensor
