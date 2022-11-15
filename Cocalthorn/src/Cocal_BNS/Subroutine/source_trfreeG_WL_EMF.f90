subroutine source_trfreeG_WL_EMF(souten)
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf, ntgeq
  use coordinate_grav_r, only : hrg
  use interface_sourceterm_trfreeG_WL
  use interface_sourceterm_trfreeG_WL_EMF
  use interface_sourceterm_trfreeG_WL_SEM
  use interface_interpo_fl2gr_midpoint
  use interface_correct_matter_source_midpoint
  use make_array_3d
  use make_array_4d
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8), pointer :: soutenf(:,:,:,:), souf(:,:,:), soug(:,:,:)
  real(8), pointer :: sou1(:,:,:,:), sou2(:,:,:,:), sou3(:,:,:,:)
  integer :: ia, irg, itg, ipg
!
  call alloc_array4d(sou1,0,nrg,0,ntg,0,npg,1,6)
  call alloc_array4d(sou2,0,nrg,0,ntg,0,npg,1,6)
  call alloc_array4d(sou3,0,nrg,0,ntg,0,npg,1,6)
  call alloc_array4d(soutenf,0,nrf,0,ntf,0,npf,1,6)
  call alloc_array3d(souf,0,nrf,0,ntf,0,npf)
  call alloc_array3d(soug,0,nrg,0,ntg,0,npg)
!
  call sourceterm_trfreeG_WL(sou1)
  call sourceterm_trfreeG_WL_EMF(sou2)
!bad routine
  call sourceterm_trfreeG_WL_SEM(soutenf)
!bad routine
  do ia = 1, 6
    souf(0:nrf,0:ntf,0:npf) = soutenf(0:nrf,0:ntf,0:npf,ia)
    call interpo_fl2gr_midpoint(souf,soug)
    call correct_matter_source_midpoint(soug)
    sou3(0:nrg,0:ntg,0:npg,ia) = soug(0:nrg,0:ntg,0:npg)
  end do
  souten(0:nrg,0:ntg,0:npg,1:6) = sou1(0:nrg,0:ntg,0:npg,1:6) &
  &                             + sou2(0:nrg,0:ntg,0:npg,1:6) &
  &                             + sou3(0:nrg,0:ntg,0:npg,1:6)
!
!!      itg = ntgeq; ipg = npgxzp
!      itg = 73; ipg = 16
!      open(15,file='test_vec_sou',status='unknown')
!        do irg = 1, nrg
!          write(15,'(1p,12e20.12)')  hrg(irg),  sou1(irg,itg,ipg,1) &
!              &                              , sou2(irg,itg,ipg,1) &
!              &                              , sou3(irg,itg,ipg,1) &
!              &                              , sou1(irg,itg,ipg,2) &
!              &                              , sou2(irg,itg,ipg,2) &
!              &                              , sou3(irg,itg,ipg,2) &
!              &                              , sou1(irg,itg,ipg,4) &
!              &                              , sou2(irg,itg,ipg,4) &
!              &                              , sou3(irg,itg,ipg,4)
!        end do
!      close(15)
!!
  deallocate(sou1)
  deallocate(sou2)
  deallocate(sou3)
  deallocate(soutenf)
  deallocate(souf)
  deallocate(soug)
!
end subroutine source_trfreeG_WL_EMF
