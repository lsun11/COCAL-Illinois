module grid_points_binary_excision_mpt
  use phys_constant, only : long, nnmpt
  implicit none
  integer, pointer :: irg_exin_(:,:,:),  itg_exin_(:,:,:),  ipg_exin_(:,:,:)
  integer, pointer :: irg_exout_(:,:,:), itg_exout_(:,:,:), ipg_exout_(:,:,:)
  real(long), pointer :: rg_exin_(:,:,:), thg_exin_(:,:,:), phig_exin_(:,:,:)
  real(long), pointer :: rg_exout_(:,:,:),thg_exout_(:,:,:),phig_exout_(:,:,:)
  integer, pointer :: ihrg_exin_(:,:,:), ihtg_exin_(:,:,:), ihpg_exin_(:,:,:)
  integer, pointer :: ihrg_exout_(:,:,:),ihtg_exout_(:,:,:),ihpg_exout_(:,:,:)
  real(long), pointer :: hrg_exin_(:,:,:),hthg_exin_(:,:,:),hphig_exin_(:,:,:)
  real(long), pointer :: hrg_exout_(:,:,:),hthg_exout_(:,:,:), &
  &                      hphig_exout_(:,:,:)
!
  real(long), pointer :: rb_(:,:,:,:), thb_(:,:,:,:), phib_(:,:,:,:)
  real(long), pointer :: hrb_(:,:,:,:), hthb_(:,:,:,:), hphib_(:,:,:,:)
  integer :: ntg_exin_min_(nnmpt), ntg_exout_max_(nnmpt), &
  &          npg_exin_max_(nnmpt), npg_exout_min_(nnmpt)
end module grid_points_binary_excision_mpt
