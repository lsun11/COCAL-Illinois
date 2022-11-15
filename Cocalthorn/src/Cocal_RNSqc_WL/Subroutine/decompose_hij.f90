subroutine decompose_hij
  use phys_constant, only : pi, long
  use def_bh_parameter
  use def_kerr_schild
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
          &                  qxxd, qxyd, qxzd, qyyd, qyzd, qzzd
  implicit none
!
  hxxd(0:nrg,0:ntg,0:npg) = hxxd_ks(0:nrg,0:ntg,0:npg) + qxxd(0:nrg,0:ntg,0:npg)
  hxyd(0:nrg,0:ntg,0:npg) = hxyd_ks(0:nrg,0:ntg,0:npg) + qxyd(0:nrg,0:ntg,0:npg)
  hxzd(0:nrg,0:ntg,0:npg) = hxzd_ks(0:nrg,0:ntg,0:npg) + qxzd(0:nrg,0:ntg,0:npg)
  hyyd(0:nrg,0:ntg,0:npg) = hyyd_ks(0:nrg,0:ntg,0:npg) + qyyd(0:nrg,0:ntg,0:npg)
  hyzd(0:nrg,0:ntg,0:npg) = hyzd_ks(0:nrg,0:ntg,0:npg) + qyzd(0:nrg,0:ntg,0:npg)
  hzzd(0:nrg,0:ntg,0:npg) = hzzd_ks(0:nrg,0:ntg,0:npg) + qzzd(0:nrg,0:ntg,0:npg)
!
end subroutine decompose_hij
