subroutine copy_grid_points_binary_excision_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_binary_excision
  use grid_points_binary_excision_mpt
  use copy_array_2dto3d_mpt
  use copy_array_3dto4d_mpt
  use copy_int_array_2dto3d_mpt
  use copy_int_array_static_0dto1d_mpt
  implicit none
  integer :: impt
!
  call copy_int_array2dto3d_mpt(impt, irg_exin, irg_exin_, 0, ntg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, irg_exout, irg_exout_, 0, ntg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, itg_exin, itg_exin_, 0, nrg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, itg_exout, itg_exout_, 0, nrg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, ipg_exin, ipg_exin_, 0, nrg, 0, ntg)
  call copy_int_array2dto3d_mpt(impt, ipg_exout, ipg_exout_, 0, nrg, 0, ntg)
!
  call copy_array2dto3d_mpt(impt, rg_exin, rg_exin_, 0, ntg, 0, npg)
  call copy_array2dto3d_mpt(impt, rg_exout, rg_exout_, 0, ntg, 0, npg)
  call copy_array2dto3d_mpt(impt, thg_exin, thg_exin_, 0, nrg, 0, npg)
  call copy_array2dto3d_mpt(impt, thg_exout, thg_exout_, 0, nrg, 0, npg)
  call copy_array2dto3d_mpt(impt, phig_exin, phig_exin_, 0, nrg, 0, ntg)
  call copy_array2dto3d_mpt(impt, phig_exout, phig_exout_, 0, nrg, 0, ntg)
!
  call copy_int_array2dto3d_mpt(impt, ihrg_exin, ihrg_exin_, 0, ntg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, ihrg_exout, ihrg_exout_, 0, ntg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, ihtg_exin, ihtg_exin_, 0, nrg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, ihtg_exout, ihtg_exout_, 0, nrg, 0, npg)
  call copy_int_array2dto3d_mpt(impt, ihpg_exin, ihpg_exin_, 0, nrg, 0, ntg)
  call copy_int_array2dto3d_mpt(impt, ihpg_exout, ihpg_exout_, 0, nrg, 0, ntg)
!
  call copy_array2dto3d_mpt(impt, hrg_exin, hrg_exin_, 0, ntg, 0, npg)
  call copy_array2dto3d_mpt(impt, hrg_exout, hrg_exout_,0, ntg, 0, npg)
  call copy_array2dto3d_mpt(impt, hthg_exin, hthg_exin_, 0, nrg, 0, npg)
  call copy_array2dto3d_mpt(impt, hthg_exout, hthg_exout_, 0, nrg, 0, npg)
  call copy_array2dto3d_mpt(impt, hphig_exin, hphig_exin_, 0, nrg, 0, ntg)
  call copy_array2dto3d_mpt(impt, hphig_exout, hphig_exout_, 0, nrg, 0, ntg)
!
  call copy_array3dto4d_mpt(impt, rb, rb_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, thb, thb_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, phib, phib_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, hrb, hrb_, 1, nrg, 1, ntg, 1, npg)
  call copy_array3dto4d_mpt(impt, hthb, hthb_, 1, nrg, 1, ntg, 1, npg)
  call copy_array3dto4d_mpt(impt, hphib, hphib_, 1, nrg, 1, ntg, 1, npg)
!
  call copy_int_arraystatic_0dto1d_mpt(impt, ntg_exin_min, ntg_exin_min_)
  call copy_int_arraystatic_0dto1d_mpt(impt, ntg_exout_max, ntg_exout_max_)
  call copy_int_arraystatic_0dto1d_mpt(impt, npg_exin_max, npg_exin_max_)
  call copy_int_arraystatic_0dto1d_mpt(impt, npg_exout_min, npg_exout_min_)
!
end subroutine copy_grid_points_binary_excision_to_mpt
