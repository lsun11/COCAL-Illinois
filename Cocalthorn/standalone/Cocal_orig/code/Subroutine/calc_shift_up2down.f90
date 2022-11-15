subroutine calc_shift_up2down
  use def_metric, only : bvxu, bvyu, bvzu, bvxd, bvyd, bvzd
  use interface_index_vec_up2down
  implicit none
  call index_vec_up2down(bvxd,bvyd,bvzd,bvxu,bvyu,bvzu)
end subroutine calc_shift_up2down
