subroutine copy_def_quantities_bh_from_mpt(impt)
  use def_quantities_bh
  use def_quantities_bh_mpt
  implicit none
  integer :: i, impt
!
  i=0
  i=i+1; AHarea = def_quantities_bh_(i,impt)
  i=i+1; AHmass = def_quantities_bh_(i,impt)
  i=i+1; AHspin = def_quantities_bh_(i,impt)
!
end subroutine copy_def_quantities_bh_from_mpt
