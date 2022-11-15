module def_SEM_tensor_EMF
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  rhoH_EMF(:,:,:)
  real(long), pointer  ::   jmd_EMF(:,:,:,:),     jmu_EMF(:,:,:,:)
  real(long), pointer  :: smijd_EMF(:,:,:,:,:), smiju_EMF(:,:,:,:,:)
  real(long), pointer  ::  trsm_EMF(:,:,:)
  real(long), pointer  ::  rhoH_EMF_grid(:,:,:)
  real(long), pointer  ::   jmd_EMF_grid(:,:,:,:),     jmu_EMF_grid(:,:,:,:)
  real(long), pointer  :: smijd_EMF_grid(:,:,:,:,:), smiju_EMF_grid(:,:,:,:,:)
  real(long), pointer  ::  trsm_EMF_grid(:,:,:)
end module def_SEM_tensor_EMF
