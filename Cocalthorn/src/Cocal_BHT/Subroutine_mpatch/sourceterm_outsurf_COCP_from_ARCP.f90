subroutine sourceterm_outsurf_COCP_from_ARCP(impt,fnchar,char,sou_out,dsou_out)
  use phys_constant, only :  long, nmpt
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use interface_sourceterm_outsurf_interpo_from_asympto_parity_mpt
  use def_metric_mpt, only : psi_, alps_, bvxd_, bvyd_, bvzd_
  use copy_array_4dto3d_mpt
  use make_array_3d
  implicit none
  real(long), pointer :: fnc_itp(:,:,:)
  real(long), pointer :: sou_out(:,:), dsou_out(:,:)
  integer :: impt
  character(len=4) :: fnchar
  character(len=2) :: char
!
  call copy_grid_parameter_interpo_from_mpt(nmpt)
  call copy_coordinate_grav_extended_interpo_from_mpt(nmpt)
  call alloc_array3d(fnc_itp,0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'psi ') call copy_array4dto3d_mpt(nmpt,psi_,fnc_itp, &
  &                                   0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'alps') call copy_array4dto3d_mpt(nmpt,alps_,fnc_itp, &
  &                                    0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'bvxd') call copy_array4dto3d_mpt(nmpt,bvxd_,fnc_itp, &
  &                                    0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'bvyd') call copy_array4dto3d_mpt(nmpt,bvyd_,fnc_itp, &
  &                                    0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'bvzd') call copy_array4dto3d_mpt(nmpt,bvzd_,fnc_itp, &
  &                                    0,nrg_itp,0,ntg_itp,0,npg_itp)
  call copy_grid_points_binary_in_asympto_from_mpt(impt)
  call sourceterm_outsurf_interpo_from_asympto_parity_mpt &
  &                              (impt,nmpt,fnc_itp,sou_out,dsou_out,char)
!
  deallocate(fnc_itp)
!
end subroutine sourceterm_outsurf_COCP_from_ARCP
