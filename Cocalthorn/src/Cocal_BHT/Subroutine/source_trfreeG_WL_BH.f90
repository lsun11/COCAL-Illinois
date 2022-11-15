subroutine source_trfreeG_WL_BH(souten)
  use grid_parameter, only : nrg, ntg, npg
  use interface_sourceterm_trfreeG_WL
  use make_array_4d
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8), pointer :: sou2(:,:,:,:)
  integer :: ia
!
  call alloc_array4d(sou2,0,nrg,0,ntg,0,npg,1,6)
!
  call sourceterm_trfreeG_WL(sou2)
  souten(0:nrg,0:ntg,0:npg,1:6) = sou2(0:nrg,0:ntg,0:npg,1:6)
!
  deallocate(sou2)
!
end subroutine source_trfreeG_WL_BH
