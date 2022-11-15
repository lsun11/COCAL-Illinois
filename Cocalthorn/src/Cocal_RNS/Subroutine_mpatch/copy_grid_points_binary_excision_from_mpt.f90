subroutine copy_grid_points_binary_excision_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_binary_excision
  use grid_points_binary_excision_mpt
  use copy_array_3dto2d_mpt
  use copy_array_4dto3d_mpt
  use copy_int_array_3dto2d_mpt
  use copy_int_array_static_1dto0d_mpt
  implicit none
  integer :: impt
!
  call copy_int_array3dto2d_mpt(impt, irg_exin_, irg_exin, 0, ntg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, irg_exout_, irg_exout, 0, ntg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, itg_exin_, itg_exin, 0, nrg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, itg_exout_, itg_exout, 0, nrg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, ipg_exin_, ipg_exin, 0, nrg, 0, ntg)
  call copy_int_array3dto2d_mpt(impt, ipg_exout_, ipg_exout, 0, nrg, 0, ntg)
!
  call copy_array3dto2d_mpt(impt, rg_exin_, rg_exin, 0, ntg, 0, npg)
  call copy_array3dto2d_mpt(impt, rg_exout_, rg_exout, 0, ntg, 0, npg)
  call copy_array3dto2d_mpt(impt, thg_exin_, thg_exin, 0, nrg, 0, npg)
  call copy_array3dto2d_mpt(impt, thg_exout_, thg_exout, 0, nrg, 0, npg)
  call copy_array3dto2d_mpt(impt, phig_exin_, phig_exin, 0, nrg, 0, ntg)
  call copy_array3dto2d_mpt(impt, phig_exout_, phig_exout, 0, nrg, 0, ntg)
!
  call copy_int_array3dto2d_mpt(impt, ihrg_exin_, ihrg_exin, 0, ntg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, ihrg_exout_, ihrg_exout, 0, ntg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, ihtg_exin_, ihtg_exin, 0, nrg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, ihtg_exout_, ihtg_exout, 0, nrg, 0, npg)
  call copy_int_array3dto2d_mpt(impt, ihpg_exin_, ihpg_exin, 0, nrg, 0, ntg)
  call copy_int_array3dto2d_mpt(impt, ihpg_exout_, ihpg_exout, 0, nrg, 0, ntg)
!
  call copy_array3dto2d_mpt(impt, hrg_exin_, hrg_exin, 0, ntg, 0, npg)
  call copy_array3dto2d_mpt(impt, hrg_exout_, hrg_exout,0, ntg, 0, npg)
  call copy_array3dto2d_mpt(impt, hthg_exin_, hthg_exin, 0, nrg, 0, npg)
  call copy_array3dto2d_mpt(impt, hthg_exout_, hthg_exout, 0, nrg, 0, npg)
  call copy_array3dto2d_mpt(impt, hphig_exin_, hphig_exin, 0, nrg, 0, ntg)
  call copy_array3dto2d_mpt(impt, hphig_exout_, hphig_exout, 0, nrg, 0, ntg)
!
  call copy_array4dto3d_mpt(impt, rb_, rb, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, thb_, thb, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, phib_, phib, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, hrb_, hrb, 1, nrg, 1, ntg, 1, npg)
  call copy_array4dto3d_mpt(impt, hthb_, hthb, 1, nrg, 1, ntg, 1, npg)
  call copy_array4dto3d_mpt(impt, hphib_, hphib, 1, nrg, 1, ntg, 1, npg)
!
  call copy_int_arraystatic_1dto0d_mpt(impt, ntg_exin_min_, ntg_exin_min)
  call copy_int_arraystatic_1dto0d_mpt(impt, ntg_exout_max_, ntg_exout_max)
  call copy_int_arraystatic_1dto0d_mpt(impt, npg_exin_max_, npg_exin_max)
  call copy_int_arraystatic_1dto0d_mpt(impt, npg_exout_min_, npg_exout_min)
!
end subroutine copy_grid_points_binary_excision_from_mpt
