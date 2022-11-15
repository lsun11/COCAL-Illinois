subroutine reset_metric_CF
  use grid_parameter, only : nrg, ntg, npg, ntgxy, npgxzm
  use def_metric, only : bvxd, bvzd
  use def_bh_parameter, only : bh_bctype, bh_sptype
  implicit none
!
  if (bh_bctype.eq.'TU') then
    if (bh_sptype.eq.'Zp'.or.bh_sptype.eq.'Zm') then
! Clear the value of bvxd on xz-plane
      bvxd(0:nrg,0:ntg,0) = 0.0d0
      bvxd(0:nrg,0:ntg,npgxzm) = 0.0d0
      bvxd(0:nrg,0:ntg,npg) = 0.0d0
      bvxd(0:nrg,  0,0:npg) = 0.0d0
      bvxd(0:nrg,ntg,0:npg) = 0.0d0
! Clear the value of bvzd on xz-plane
      bvzd(0:nrg,0:ntg,0) = 0.0d0
      bvzd(0:nrg,0:ntg,npgxzm) = 0.0d0
      bvzd(0:nrg,0:ntg,npg) = 0.0d0
      bvzd(0:nrg,  0,0:npg) = 0.0d0
      bvzd(0:nrg,ntg,0:npg) = 0.0d0
    end if
  end if
!
  if (bh_sptype.eq.'Zp'.or.bh_sptype.eq.'Zm') then
! Clear the value of bvzd on xy-plane
    bvzd(0:nrg,ntgxy,0:npg) = 0.0d0
  end if
!
end subroutine reset_metric_CF
