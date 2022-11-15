subroutine interpolation_fillup_binary_COCP(impt,fnchar,char,fnc)
  use phys_constant, only :  long
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use interface_interpolation_fillup_binary_parity_mpt
  use def_metric_mpt, only : psi_, alph_, bvxd_, bvyd_, bvzd_
  use copy_array_4dto3d_mpt
  use make_array_3d
 implicit none
  real(long), pointer :: fnc_itp(:,:,:), fnc(:,:,:)
  integer :: impt, impt_ex, irg
  character(len=4) :: fnchar
  character(len=2) :: char
!
  if (impt.eq.1) impt_ex = 2
  if (impt.eq.2) impt_ex = 1
  call copy_grid_parameter_interpo_from_mpt(impt_ex)
  call copy_grid_parameter_binary_excision_interpo_from_mpt(impt_ex)
  call copy_coordinate_grav_extended_interpo_from_mpt(impt_ex)
  call alloc_array3d(fnc_itp,0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'psi ') call copy_array4dto3d_mpt(impt_ex,psi_,fnc_itp, &
  &                                      0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'alph') call copy_array4dto3d_mpt(impt_ex,alph_,fnc_itp, &
  &                                       0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'bvxd') call copy_array4dto3d_mpt(impt_ex,bvxd_,fnc_itp, &
  &                                       0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'bvyd') call copy_array4dto3d_mpt(impt_ex,bvyd_,fnc_itp, &
  &                                       0,nrg_itp,0,ntg_itp,0,npg_itp)
  if (fnchar.eq.'bvzd') call copy_array4dto3d_mpt(impt_ex,bvzd_,fnc_itp, &
  &                                       0,nrg_itp,0,ntg_itp,0,npg_itp)
  call interpolation_fillup_binary_parity_mpt(fnc,fnc_itp,char)
!
  deallocate(fnc_itp)
!
end subroutine interpolation_fillup_binary_COCP
