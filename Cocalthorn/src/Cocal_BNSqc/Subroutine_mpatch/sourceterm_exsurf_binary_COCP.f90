subroutine sourceterm_exsurf_binary_COCP(impt,fnchar,char,sou_ex,dsou_ex)
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use interface_sourceterm_exsurf_binary_parity
  use def_metric, only : psi, alps, bvxd, bvyd, bvzd
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
  call copy_def_metric_from_mpt(impt_ex)
  call alloc_array3d(fnc,0,nrg,0,ntg,0,npg)
  if (fnchar.eq.'psi ') fnc(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)
  if (fnchar.eq.'alps') fnc(0:nrg,0:ntg,0:npg) = alps(0:nrg,0:ntg,0:npg)
  if (fnchar.eq.'bvxd') fnc(0:nrg,0:ntg,0:npg) = bvxd(0:nrg,0:ntg,0:npg)
  if (fnchar.eq.'bvyd') fnc(0:nrg,0:ntg,0:npg) = bvyd(0:nrg,0:ntg,0:npg)
  if (fnchar.eq.'bvzd') fnc(0:nrg,0:ntg,0:npg) = bvzd(0:nrg,0:ntg,0:npg)
  call sourceterm_exsurf_binary_parity(fnc,sou_ex,dsou_ex,char)
  call copy_from_mpatch_exsurf_binary_COCP(impt)
  call copy_def_metric_from_mpt(impt)
!
  deallocate(fnc)
!
end subroutine sourceterm_exsurf_binary_COCP
