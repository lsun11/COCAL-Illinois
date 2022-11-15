subroutine compute_fnc_inversion_index(pot1,pot2,index)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  implicit none
  real(long), pointer    :: pot1(:,:,:), pot2(:,:,:)
  real(long) :: index
!
  pot2(0:nrg,0:ntg,0:npg) = pot1(0:nrg,0:ntg,0:npg)**index
end subroutine compute_fnc_inversion_index
