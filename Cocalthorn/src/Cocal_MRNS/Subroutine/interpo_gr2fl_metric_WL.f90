subroutine interpo_gr2fl_metric_WL
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd, &
  &                      bvxu, bvyu, bvzu
  use def_metric_hij
  use def_metric_on_SFC_CF
  use def_metric_on_SFC_WL
  use interface_interpo_gr2fl
  implicit none
 !
  call interpo_gr2fl(alph, alphf)
  call interpo_gr2fl(psi, psif)
  call interpo_gr2fl(bvxd, bvxdf)
  call interpo_gr2fl(bvyd, bvydf)
  call interpo_gr2fl(bvzd, bvzdf)
  call interpo_gr2fl(bvxu, bvxuf)
  call interpo_gr2fl(bvyu, bvyuf)
  call interpo_gr2fl(bvzu, bvzuf)
!
  call interpo_gr2fl(hxxd, hxxdf)
  call interpo_gr2fl(hxyd, hxydf)
  call interpo_gr2fl(hxzd, hxzdf)
  call interpo_gr2fl(hyyd, hyydf)
  call interpo_gr2fl(hyzd, hyzdf)
  call interpo_gr2fl(hzzd, hzzdf)
  call interpo_gr2fl(hxxu, hxxuf)
  call interpo_gr2fl(hxyu, hxyuf)
  call interpo_gr2fl(hxzu, hxzuf)
  call interpo_gr2fl(hyyu, hyyuf)
  call interpo_gr2fl(hyzu, hyzuf)
  call interpo_gr2fl(hzzu, hzzuf)
!
end subroutine interpo_gr2fl_metric_WL
