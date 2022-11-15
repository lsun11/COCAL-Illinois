subroutine source_trfreeG_WL_COT(souten,cobj)
  use grid_parameter, only : nrg, ntg, npg
  use def_formulation
  use interface_sourceterm_trfreeG_WL_BHT
  use interface_sourceterm_trfreeG_WL
  use interface_sourceterm_trfreeG_WL_bhex
  use make_array_4d
  implicit none
  real(8), pointer :: souten(:,:,:,:)
  real(8), pointer :: souten1(:,:,:,:), souten2(:,:,:,:)
  character(len=2), intent(in) :: cobj
!
  if (sw_disk==0) then
    if (cobj=='bh')  then
      call sourceterm_trfreeG_WL_bhex(souten)
    else
      call sourceterm_trfreeG_WL(souten)  !!!! NEED FIX
    end if
  else
    call alloc_array4d(souten1,0,nrg,0,ntg,0,npg,1,6)
    call alloc_array4d(souten2,0,nrg,0,ntg,0,npg,1,6)

    if (cobj=='bh')  then
      call sourceterm_trfreeG_WL_bhex(souten1)
      call sourceterm_trfreeG_WL_BHT(souten2)
    else
      call sourceterm_trfreeG_WL(souten1)   !!!  NEED  FIX
    end if

    souten(0:nrg,0:ntg,0:npg,1:6) = souten1(0:nrg,0:ntg,0:npg,1:6) &
    &                             + souten2(0:nrg,0:ntg,0:npg,1:6)
    deallocate(souten1)
    deallocate(souten2)
  end if

!!  call sourceterm_trfreeG_drot_SFC(soutenf)

!!  call sourceterm_trfreeG_WL(sou2)
!=>  call sourceterm_trfreeG_WL_bhex(souten)

!
end subroutine source_trfreeG_WL_COT
