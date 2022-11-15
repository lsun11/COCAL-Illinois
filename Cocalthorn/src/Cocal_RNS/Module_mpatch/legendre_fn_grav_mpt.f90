!  Associated Legendre function and factorials
!______________________________________________
module legendre_fn_grav_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer ::    plmg_(:,:,:,:),    hplmg_(:,:,:,:)
  real(long), pointer ::  dtplmg_(:,:,:,:),  hdtplmg_(:,:,:,:)
  real(long), pointer ::   yplmg_(:,:,:,:),   hyplmg_(:,:,:,:)
  real(long), pointer :: dtyplmg_(:,:,:,:), hdtyplmg_(:,:,:,:)
  real(long), pointer :: facnmg_(:,:,:), epsig_(:,:) 
end module legendre_fn_grav_mpt
