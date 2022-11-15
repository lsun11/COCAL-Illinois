subroutine copy_def_quantities_bh_to_mpt(impt)
  use def_quantities_bh
  use def_quantities_bh_mpt
  implicit none
  integer :: i, impt
!
  i=0
  i=i+1; def_quantities_bh_(i,impt) = AHarea
  i=i+1; def_quantities_bh_(i,impt) = AHmass
  i=i+1; def_quantities_bh_(i,impt) = AHspin
!
end subroutine copy_def_quantities_bh_to_mpt
