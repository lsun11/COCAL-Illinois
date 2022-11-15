subroutine update_grfield(pot,grfield,convf)
  use phys_constant,  only : long, nnrg, nntg, nnpg
  use grid_parameter, only : nrg, ntg, npg
  use make_array_3d
  implicit none
  real(long), pointer    :: pot(:,:,:)
  real(long), pointer    :: grfield(:,:,:)
  real(long), intent(in) :: convf
  real(long) :: work(0:nnrg,0:nntg,0:nnpg) = 0.0d0
!
  work(0:nrg,0:ntg,0:npg) = convf*pot(0:nrg,0:ntg,0:npg) + (1.0d0-convf)*grfield(0:nrg,0:ntg,0:npg)
  grfield(0:nrg,0:ntg,0:npg) = work(0:nrg,0:ntg,0:npg)
!
end subroutine update_grfield
