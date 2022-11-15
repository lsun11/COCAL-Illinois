subroutine update_grfield(pot,grfield,convf)
  use phys_constant,  only : long
  implicit none
  real(long), pointer    :: pot(:,:,:)
  real(long), pointer    :: grfield(:,:,:)
  real(long), intent(in) :: convf
!
  grfield = convf*pot + (1.0d0-convf)*grfield
!
end subroutine update_grfield
