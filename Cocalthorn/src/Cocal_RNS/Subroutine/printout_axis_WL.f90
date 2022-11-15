subroutine printout_axis_WL(filename,ita,ipa)
  use phys_constant, only : long
  use def_metric
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_matter
  use def_velocity_potential
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only  : nrg, ntg, npg 
  implicit none
  integer :: ita, ipa
  integer :: irg, itg, ipg,  npg_l, npg_r
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/), char_1
  character(30) :: filename
!
  npg_l = ipa+npg/2
  npg_r = ipa
  if (npg_l > npg) then
    npg_l = ipa
    npg_r = npg_l - npg
  end if
!
!  write (6,*) "npg_l, npg_r", npg_l, npg_r
!
  open(12,file=filename,status='unknown')
  do irg = nrg, 0, -1
    char_1 = ' '
    write(12,'(a1,1p,15e20.12)') char_1,  -rg(irg)   &
    &                                  ,  psi(irg,ita,npg_l) &
    &                                  , alph(irg,ita,npg_l) &
    &                                  , bvxd(irg,ita,npg_l) &
    &                                  , bvyd(irg,ita,npg_l) &
    &                                  , bvzd(irg,ita,npg_l) &
    &                                  , hxxd(irg,ita,npg_l) &
    &                                  , hxyd(irg,ita,npg_l) &
    &                                  , hxzd(irg,ita,npg_l) &
    &                                  , hyyd(irg,ita,npg_l) &
    &                                  , hyzd(irg,ita,npg_l) &
    &                                  , hzzd(irg,ita,npg_l) 
  end do
  write(12,*)  ""
  do irg = 0, nrg
    char_1 = ' '
    write(12,'(a1,1p,15e20.12)') char_1,   rg(irg)   &
    &                                  ,  psi(irg,ita,npg_r) &
    &                                  , alph(irg,ita,npg_r) &
    &                                  , bvxd(irg,ita,npg_r) &
    &                                  , bvyd(irg,ita,npg_r) &
    &                                  , bvzd(irg,ita,npg_r) &
    &                                  , hxxd(irg,ita,npg_r) &
    &                                  , hxyd(irg,ita,npg_r) &
    &                                  , hxzd(irg,ita,npg_r) &
    &                                  , hyyd(irg,ita,npg_r) &
    &                                  , hyzd(irg,ita,npg_r) &
    &                                  , hzzd(irg,ita,npg_r) 
  end do
  close(12)
!
end subroutine printout_axis_WL
