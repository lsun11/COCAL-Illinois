subroutine sourceterm_insurf_ARCP_from_COCP_pBH(impt,fnchar,parchar,sou_in, &
  &                                                                dsou_in)
  use phys_constant, only :  long, nmpt
  use grid_parameter, only : ntg, npg
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use interface_sourceterm_insurf_asympto_interpo_from_mpt
  use def_metric_pBH_mpt, only : log_wme_, log_N_
  use copy_array_4dto3d_mpt
  use make_array_2d
  use make_array_3d
  implicit none
  real(long), pointer :: fnc_itp(:,:,:)
  real(long), pointer :: sou_in(:,:), dsou_in(:,:)
  real(long), pointer :: sou_work(:,:), dsou_work(:,:)
  real(long) :: pari
  integer :: impt, impt_bin
  character(len=4) :: fnchar
  character(len=2) :: parchar
!
  call alloc_array2d(sou_work,0,ntg,0,npg)
  call alloc_array2d(dsou_work,0,ntg,0,npg)
!
  if (parchar.eq.'ev') pari = + 1.0d0
  if (parchar.eq.'od') pari = - 1.0d0
  do impt_bin = 1, 2
    call copy_grid_parameter_interpo_from_mpt(impt_bin)
    call copy_coordinate_grav_extended_interpo_from_mpt(impt_bin)
    call alloc_array3d(fnc_itp,0,nrg_itp,0,ntg_itp,0,npg_itp)
    if (fnchar.eq.'logw') call copy_array4dto3d_mpt(impt_bin,log_wme_, &
    &                                 fnc_itp,0,nrg_itp,0,ntg_itp,0,npg_itp)
    if (fnchar.eq.'logN') call copy_array4dto3d_mpt(impt_bin,log_N_, &
    &                                 fnc_itp,0,nrg_itp,0,ntg_itp,0,npg_itp)
    call copy_grid_points_asymptotic_patch_from_mpt(impt_bin)
    call sourceterm_insurf_asympto_interpo_from_mpt &
    &                           (impt_bin,impt,fnc_itp,sou_in,dsou_in)
    deallocate(fnc_itp)
    if (impt_bin.eq.1) then
       sou_work(0:ntg,0:npg) =  sou_in(0:ntg,0:npg)
      dsou_work(0:ntg,0:npg) = dsou_in(0:ntg,0:npg)
    end if
    if (impt_bin.eq.2) then
       sou_in(0:ntg,0:npg) = &
      &   0.5d0*( sou_work(0:ntg,0:npg) + pari* sou_in(0:ntg,0:npg))
      dsou_in(0:ntg,0:npg) = &
      &   0.5d0*(dsou_work(0:ntg,0:npg) + pari*dsou_in(0:ntg,0:npg))
    end if
  end do
!
  deallocate(sou_work)
  deallocate(dsou_work)
!
end subroutine sourceterm_insurf_ARCP_from_COCP_pBH
