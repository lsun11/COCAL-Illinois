subroutine IO_output_1D_xyz(iter,fch)
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_matter, only : emdg, omeg, utg
  use def_CTT_decomposition, only : wvxd, wvyd, wvzd
  use interface_IO_output_1D_general
  use make_array_3d
  implicit none
  real(long), pointer :: aps(:,:,:)
  integer :: its, ips
  character(30) :: char3,char2,char1
  integer, intent(in) :: iter
  character(len=2), intent(in) :: fch
!
  call alloc_array3d(aps,0,nrg,0,ntg,0,npg)

  its = 6
  ips = 3

  write(char1, '(i5)') iter
  char2 = adjustl(char1)

  char3 = 'psi_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', psi, -1, ntg/2,0)
  char3 = 'psi_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', psi, -1, ntg/2,npg/4)
  char3 = 'psi_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', psi, -1, 0,0)
!  char3 = 'psi_ar' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', psi, -1, its,ips)

  aps(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)*alph(0:nrg,0:ntg,0:npg)
  char3 = 'alps_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', aps, -1, ntg/2,0)
  char3 = 'alps_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', aps, -1, ntg/2,npg/4)
  char3 = 'alps_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', aps, -1, 0,0)
!  char3 = 'alps_ar' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', aps, -1, its,ips)

  char3 = 'bvxd_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvxd, -1, ntg/2,0)
  char3 = 'bvxd_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvxd, -1, ntg/2,npg/4)
  char3 = 'bvxd_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvxd, -1, 0,0)
!  char3 = 'bvxd_ar' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', bvxd, -1, its,ips)

  char3 = 'bvyd_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvyd, -1, ntg/2,0)
  char3 = 'bvyd_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvyd, -1, ntg/2,npg/4)
  char3 = 'bvyd_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvyd, -1, 0,0)
!  char3 = 'bvyd_ar' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', bvyd, -1, its,ips)

  char3 = 'bvzd_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvzd, -1, ntg/2,0)
  char3 = 'bvzd_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvzd, -1, ntg/2,npg/4)
  char3 = 'bvzd_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', bvzd, -1, 0,0)
!  char3 = 'bvzd_ar' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', bvzd, -1, its,ips)

  char3 = 'emdg_xp' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', emdg, -1, ntg/2,0)
  char3 = 'emdg_xn' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', emdg, -1, ntg/2,npg/2)
  char3 = 'emdg_yp' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', emdg, -1, ntg/2,npg/4)
  char3 = 'emdg_yn' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', emdg, -1, ntg/2,3*npg/4)

  char3 = 'omeg_xp' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', omeg, -1, ntg/2,0)
  char3 = 'omeg_xn' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', omeg, -1, ntg/2,npg/2)
  char3 = 'omeg_yp' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', omeg, -1, ntg/2,npg/4)
  char3 = 'omeg_yn' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', omeg, -1, ntg/2,3*npg/4)

  char3 = 'utg_xp' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', utg, -1, ntg/2,0)
  char3 = 'utg_xn' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', utg, -1, ntg/2,npg/2)
  char3 = 'utg_yp' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', utg, -1, ntg/2,npg/4)
  char3 = 'utg_yn' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', utg, -1, ntg/2,3*npg/4)

  if (fch=='WL') then
    char3 = 'hxxd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxxd, -1, ntg/2,0)
    char3 = 'hxxd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxxd, -1, ntg/2,npg/4)
    char3 = 'hxxd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxxd, -1, 0,0)
!    char3 = 'hxxd_ar' // trim(char2) // '.txt'
!    call IO_output_1D_general(char3, 'g', 'g', hxxd, -1, its,ips)
  
    char3 = 'hxyd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxyd, -1, ntg/2,0)
    char3 = 'hxyd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxyd, -1, ntg/2,npg/4)
    char3 = 'hxyd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxyd, -1, 0,0)
!    char3 = 'hxyd_ar' // trim(char2) // '.txt'
!    call IO_output_1D_general(char3, 'g', 'g', hxyd, -1, its,ips)
  
    char3 = 'hxzd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxzd, -1, ntg/2,0)
    char3 = 'hxzd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxzd, -1, ntg/2,npg/4)
    char3 = 'hxzd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hxzd, -1, 0,0)
!    char3 = 'hxzd_ar' // trim(char2) // '.txt'
!    call IO_output_1D_general(char3, 'g', 'g', hxzd, -1, its,ips)
  
    char3 = 'hyyd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hyyd, -1, ntg/2,0)
    char3 = 'hyyd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hyyd, -1, ntg/2,npg/4)
    char3 = 'hyyd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hyyd, -1, 0,0)
!    char3 = 'hyyd_ar' // trim(char2) // '.txt'
!    call IO_output_1D_general(char3, 'g', 'g', hyyd, -1, its,ips)
  
    char3 = 'hyzd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hyzd, -1, ntg/2,0)
    char3 = 'hyzd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hyzd, -1, ntg/2,npg/4)
    char3 = 'hyzd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hyzd, -1, 0,0)
!    char3 = 'hyzd_ar' // trim(char2) // '.txt'
!    call IO_output_1D_general(char3, 'g', 'g', hyzd, -1, its,ips)
  
    char3 = 'hzzd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hzzd, -1, ntg/2,0)
    char3 = 'hzzd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hzzd, -1, ntg/2,npg/4)
    char3 = 'hzzd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', hzzd, -1, 0,0)
!    char3 = 'hzzd_ar' // trim(char2) // '.txt'
!    call IO_output_1D_general(char3, 'g', 'g', hzzd, -1, its,ips)

    char3 = 'wvxd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvxd, -1, ntg/2,0)
    char3 = 'wvxd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvxd, -1, ntg/2,npg/4)
    char3 = 'wvxd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvxd, -1, 0,0)

    char3 = 'wvyd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvyd, -1, ntg/2,0)
    char3 = 'wvyd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvyd, -1, ntg/2,npg/4)
    char3 = 'wvyd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvyd, -1, 0,0)

    char3 = 'wvzd_xa' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvzd, -1, ntg/2,0)
    char3 = 'wvzd_ya' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvzd, -1, ntg/2,npg/4)
    char3 = 'wvzd_za' // trim(char2) // '.txt'
    call IO_output_1D_general(char3, 'g', 'g', wvzd, -1, 0,0)
  end if  
!
  deallocate(aps)
!
end subroutine IO_output_1D_xyz
