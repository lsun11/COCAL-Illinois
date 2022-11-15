subroutine copy_coordinate_grav_extended_interpo_from_mpt(impt)
  use phys_constant, only : long, nnrg, nntg, nnpg
!  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use coordinate_grav_extended_interpo
  use coordinate_grav_extended_mpt
  use copy_array_static_2dto1d_mpt
  use copy_int_array_static_2dto1d_mpt
  use copy_int_array_static_3dto2d_mpt
  implicit none
  integer :: impt
!
  call copy_arraystatic_2dto1d_mpt(impt,rgex_,rgex_itp,-2,nnrg+2)
  call copy_arraystatic_2dto1d_mpt(impt,thgex_,thgex_itp,-2,nntg+2)
  call copy_arraystatic_2dto1d_mpt(impt,phigex_,phigex_itp,-2,nnpg+2)
  call copy_int_arraystatic_2dto1d_mpt(impt,irgex_r_,irgex_r_itp,-2,nnrg+2)
  call copy_int_arraystatic_3dto2d_mpt(impt,itgex_r_,itgex_r_itp, 0,nntg,-2,nnrg+2)
  call copy_int_arraystatic_3dto2d_mpt(impt,ipgex_r_,ipgex_r_itp, 0,nnpg,-2,nnrg+2)
  call copy_int_arraystatic_2dto1d_mpt(impt,itgex_th_,itgex_th_itp,-2,nntg+2)
  call copy_int_arraystatic_3dto2d_mpt(impt,ipgex_th_,ipgex_th_itp, 0,nnpg,-2,nntg+2)
  call copy_int_arraystatic_2dto1d_mpt(impt,ipgex_phi_,ipgex_phi_itp,-2,nnpg+2)
!
  call copy_arraystatic_2dto1d_mpt(impt,hrgex_,hrgex_itp,-2,nnrg+2)
  call copy_arraystatic_2dto1d_mpt(impt,hthgex_,hthgex_itp,-2,nntg+2)
  call copy_arraystatic_2dto1d_mpt(impt,hphigex_,hphigex_itp,-2,nnpg+2)
  call copy_int_arraystatic_2dto1d_mpt(impt,irgex_hr_,irgex_hr_itp,-2,nnrg+2)
  call copy_int_arraystatic_3dto2d_mpt(impt,itgex_hr_,itgex_hr_itp, 1,nntg, -2,nnrg+2)
  call copy_int_arraystatic_3dto2d_mpt(impt,ipgex_hr_,ipgex_hr_itp, 1,nnpg, -2,nnrg+2)
  call copy_int_arraystatic_2dto1d_mpt(impt,itgex_hth_,itgex_hth_itp,-2,nntg+2)
  call copy_int_arraystatic_3dto2d_mpt(impt,ipgex_hth_,ipgex_hth_itp,1,nnpg,-2,nntg+2)
  call copy_int_arraystatic_2dto1d_mpt(impt,ipgex_hphi_,ipgex_hphi_itp,-2,nnpg+2)
!
end subroutine copy_coordinate_grav_extended_interpo_from_mpt
