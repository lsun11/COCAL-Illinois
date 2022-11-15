subroutine source_trfreeG_WL(souten)
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use interface_sourceterm_trfreeG_WL
  use interface_sourceterm_trfreeG_drot_SFC
  use interface_interpo_fl2gr_midpoint
  use interface_correct_matter_source_midpoint
  use make_array_3d
  use make_array_4d
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8), pointer :: soutenf(:,:,:,:), souf(:,:,:), soug(:,:,:)
  real(8), pointer :: sou1(:,:,:,:), sou2(:,:,:,:)
  integer :: ia
!
  call alloc_array4d(soutenf,0,nrf,0,ntf,0,npf,1,6)
  call alloc_array3d(souf,0,nrf,0,ntf,0,npf)
  call alloc_array3d(soug,0,nrg,0,ntg,0,npg)
  call alloc_array4d(sou1,0,nrg,0,ntg,0,npg,1,6)
  call alloc_array4d(sou2,0,nrg,0,ntg,0,npg,1,6)
!
  call sourceterm_trfreeG_drot_SFC(soutenf)
  do ia = 1, 6
    souf(0:nrf,0:ntf,0:npf) = soutenf(0:nrf,0:ntf,0:npf,ia)
    call interpo_fl2gr_midpoint(souf,soug)
    call correct_matter_source_midpoint(soug)
    sou1(0:nrg,0:ntg,0:npg,ia) = soug(0:nrg,0:ntg,0:npg)
  end do
  call sourceterm_trfreeG_WL(sou2)
  souten(0:nrg,0:ntg,0:npg,1:6) = sou1(0:nrg,0:ntg,0:npg,1:6) &
  &                             + sou2(0:nrg,0:ntg,0:npg,1:6)
!
  deallocate(soutenf)
  deallocate(souf)
  deallocate(soug)
  deallocate(sou1)
  deallocate(sou2)
!
end subroutine source_trfreeG_WL
