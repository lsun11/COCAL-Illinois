subroutine interpo_gr2fl_metric_CF
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd, &
  &                      bvxu, bvyu, bvzu
  use def_metric_on_SFC_CF
  use interface_interpo_gr2fl
  implicit none
 !
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
  call interpo_gr2fl(bvxd, bvxdf)
  call interpo_gr2fl(bvyd, bvydf)
  call interpo_gr2fl(bvzd, bvzdf)
!  call interpo_gr2fl(bvxu, bvxuf)
!  call interpo_gr2fl(bvyu, bvyuf)
!  call interpo_gr2fl(bvzu, bvzuf)
!
end subroutine interpo_gr2fl_metric_CF
