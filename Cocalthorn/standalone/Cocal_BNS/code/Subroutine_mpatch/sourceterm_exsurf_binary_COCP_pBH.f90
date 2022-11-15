subroutine sourceterm_exsurf_binary_COCP_pBH(impt,fnchar,char,sou_ex,dsou_ex)
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use interface_sourceterm_exsurf_binary_parity
  use def_metric_pBH, only : log_wme, log_N
  use make_array_3d
  implicit none
  real(long), pointer :: fnc(:,:,:), sou_ex(:,:), dsou_ex(:,:)
  integer :: impt, impt_ex
  character(len=4) :: fnchar
  character(len=2) :: char
!
  if (impt.eq.1) impt_ex = 2
  if (impt.eq.2) impt_ex = 1
  call copy_from_mpatch_exsurf_binary_COCP(impt_ex)
  call copy_def_metric_pBH_from_mpt(impt_ex)
  call alloc_array3d(fnc,0,nrg,0,ntg,0,npg)
  if (fnchar.eq.'logw') fnc(0:nrg,0:ntg,0:npg) = log_wme(0:nrg,0:ntg,0:npg)
  if (fnchar.eq.'logN') fnc(0:nrg,0:ntg,0:npg) = log_N(0:nrg,0:ntg,0:npg)
  call sourceterm_exsurf_binary_parity(fnc,sou_ex,dsou_ex,char)
  call copy_from_mpatch_exsurf_binary_COCP(impt)
  call copy_def_metric_pBH_from_mpt(impt)
!
  deallocate(fnc)
!
end subroutine sourceterm_exsurf_binary_COCP_pBH
