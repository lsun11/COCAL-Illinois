subroutine calc_shift_down2up
  use def_metric, only : bvxu, bvyu, bvzu, bvxd, bvyd, bvzd
  use interface_index_vec_down2up
  implicit none
  call index_vec_down2up(bvxu,bvyu,bvzu,bvxd,bvyd,bvzd)
end subroutine calc_shift_down2up
