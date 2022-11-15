subroutine compute_fnc_division(pot1,pot2,pot3)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  implicit none
  real(long), pointer    :: pot1(:,:,:), pot2(:,:,:), pot3(:,:,:)
!
  pot3(0:nrg,0:ntg,0:npg) = pot1(0:nrg,0:ntg,0:npg)/pot2(0:nrg,0:ntg,0:npg)
end subroutine compute_fnc_division
