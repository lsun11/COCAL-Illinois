subroutine IO_output_1D_exact_kerr_schild
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg
  use def_kerr_schild
  use interface_IO_output_1D_general
  use make_array_3d
  implicit none
  real(long), pointer :: aps(:,:,:)
  integer :: its, ips
  character(30) :: char3
!
  call alloc_array3d(aps,0,nrg,0,ntg,0,npg)

  its = 6
  ips = 3

  char3 = 'psi_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', psi_ks, -1, ntg/2,0)
  char3 = 'psi_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', psi_ks, -1, ntg/2,npg/4)
  char3 = 'psi_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', psi_ks, -1, 0,0)
!  char3 = 'psi_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', psi_ks, -1, its,ips)

  aps(0:nrg,0:ntg,0:npg) = psi_ks(0:nrg,0:ntg,0:npg)*alph_ks(0:nrg,0:ntg,0:npg)
  char3 = 'alps_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', aps, -1, ntg/2,0)
  char3 = 'alps_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', aps, -1, ntg/2,npg/4)
  char3 = 'alps_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', aps, -1, 0,0)
!  char3 = 'alps_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', aps, -1, its,ips)

  char3 = 'bvxd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvxd_ks, -1, ntg/2,0)
  char3 = 'bvxd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvxd_ks, -1, ntg/2,npg/4)
  char3 = 'bvxd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvxd_ks, -1, 0,0)
!  char3 = 'bvxd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', bvxd_ks, -1, its,ips)

  char3 = 'bvyd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvyd_ks, -1, ntg/2,0)
  char3 = 'bvyd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvyd_ks, -1, ntg/2,npg/4)
  char3 = 'bvyd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvyd_ks, -1, 0,0)
!  char3 = 'bvyd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', bvyd_ks, -1, its,ips)

  char3 = 'bvzd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvzd_ks, -1, ntg/2,0)
  char3 = 'bvzd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvzd_ks, -1, ntg/2,npg/4)
  char3 = 'bvzd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvzd_ks, -1, 0,0)
!  char3 = 'bvzd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', bvzd_ks, -1, its,ips)

  char3 = 'hxxd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxxd_ks, -1, ntg/2,0)
  char3 = 'hxxd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxxd_ks, -1, ntg/2,npg/4)
  char3 = 'hxxd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxxd_ks, -1, 0,0)
!  char3 = 'hxxd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hxxd_ks, -1, its,ips)

  char3 = 'hxyd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxyd_ks, -1, ntg/2,0)
  char3 = 'hxyd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxyd_ks, -1, ntg/2,npg/4)
  char3 = 'hxyd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxyd_ks, -1, 0,0)
!  char3 = 'hxyd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hxyd_ks, -1, its,ips)

  char3 = 'hxzd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxzd_ks, -1, ntg/2,0)
  char3 = 'hxzd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxzd_ks, -1, ntg/2,npg/4)
  char3 = 'hxzd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hxzd_ks, -1, 0,0)
!  char3 = 'hxzd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hxzd_ks, -1, its,ips)

  char3 = 'hyyd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hyyd_ks, -1, ntg/2,0)
  char3 = 'hyyd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hyyd_ks, -1, ntg/2,npg/4)
  char3 = 'hyyd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hyyd_ks, -1, 0,0)
!  char3 = 'hyyd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hyyd_ks, -1, its,ips)

  char3 = 'hyzd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hyzd_ks, -1, ntg/2,0)
  char3 = 'hyzd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hyzd_ks, -1, ntg/2,npg/4)
  char3 = 'hyzd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hyzd_ks, -1, 0,0)
!  char3 = 'hyzd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hyzd_ks, -1, its,ips)

  char3 = 'hzzd_xa_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hzzd_ks, -1, ntg/2,0)
  char3 = 'hzzd_ya_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hzzd_ks, -1, ntg/2,npg/4)
  char3 = 'hzzd_za_ex.txt'
  call IO_output_1D_general(char3, 'g', 'g', hzzd_ks, -1, 0,0)
!  char3 = 'hzzd_ar_ex.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hzzd_ks, -1, its,ips)

!
  deallocate(aps)
!
end subroutine IO_output_1D_exact_kerr_schild
