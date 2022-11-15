subroutine cleargeometry
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_cristoffel, only : cri
  use def_gamma_crist, only : gmcrix, gmcriy, gmcriz
  use def_ricci_tensor, only : rab, rabnl
  implicit none
!
  hxxd(0:nrg,0:ntg,0:npg) = 0.0d0
  hxyd(0:nrg,0:ntg,0:npg) = 0.0d0
  hxzd(0:nrg,0:ntg,0:npg) = 0.0d0
  hyyd(0:nrg,0:ntg,0:npg) = 0.0d0
  hyzd(0:nrg,0:ntg,0:npg) = 0.0d0
  hzzd(0:nrg,0:ntg,0:npg) = 0.0d0
  hxxu(0:nrg,0:ntg,0:npg) = 0.0d0
  hxyu(0:nrg,0:ntg,0:npg) = 0.0d0
  hxzu(0:nrg,0:ntg,0:npg) = 0.0d0
  hyyu(0:nrg,0:ntg,0:npg) = 0.0d0
  hyzu(0:nrg,0:ntg,0:npg) = 0.0d0
  hzzu(0:nrg,0:ntg,0:npg) = 0.0d0
  cri(0:nrg,0:ntg,0:npg,1:3,1:6) = 0.0d0
  gmcrix(0:nrg,0:ntg,0:npg) = 0.0d0
  gmcriy(0:nrg,0:ntg,0:npg) = 0.0d0
  gmcriz(0:nrg,0:ntg,0:npg) = 0.0d0
  rab(0:nrg,0:ntg,0:npg,1:6) = 0.0d0
  rabnl(0:nrg,0:ntg,0:npg,1:6) = 0.0d0
!
end subroutine cleargeometry
