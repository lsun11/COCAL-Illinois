subroutine IO_output_2D_xyz(iter,fch)
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_matter, only : emdg, omeg, utg
  use interface_IO_output_2D_general
  use make_array_3d
  implicit none
  real(long), pointer :: aps(:,:,:)
  integer :: its, ips
  character(30) :: char3,char2,char1
  integer, intent(in) :: iter
  character(len=2), intent(in) :: fch
!
  call alloc_array3d(aps,0,nrg,0,ntg,0,npg)

  write(char1, '(i5)') iter
  char2 = adjustl(char1)

  char3 = 'psi_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', psi,'xy')
  char3 = 'psi_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', psi,'xz')

  char3 = 'alph_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', alph,'xy')
  char3 = 'alph_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', alph,'xz')

  char3 = 'bvxd_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', bvxd,'xy')
  char3 = 'bvxd_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', bvxd,'xz')

  char3 = 'bvyd_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', bvyd,'xy')
  char3 = 'bvyd_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', bvyd,'xz')

  char3 = 'bvzd_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', bvzd,'xy')
  char3 = 'bvzd_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', bvzd,'xz')

  char3 = 'emdg_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', emdg,'xy')
  char3 = 'emdg_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', emdg,'xz')

  char3 = 'omeg_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', omeg,'xy')
  char3 = 'omeg_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', omeg,'xz')

  char3 = 'utg_xy' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', utg,'xy')
  char3 = 'utg_xz' // trim(char2) // '.txt'
  call IO_output_2D_general(char3, 'g', 'g', utg,'xz')

  if (fch=='WL') then
    char3 = 'hxxd_xy' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hxxd,'xy')
    char3 = 'hxxd_xz' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hxxd,'xz')
  
    char3 = 'hxyd_xy' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hxyd,'xy')
    char3 = 'hxyd_xz' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hxyd,'xz')
  
    char3 = 'hxzd_xy' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hxzd,'xy')
    char3 = 'hxzd_xz' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hxzd,'xz')
  
    char3 = 'hyyd_xy' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hyyd,'xy')
    char3 = 'hyyd_xz' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hyyd,'xz')
  
    char3 = 'hyzd_xy' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hyzd,'xy')
    char3 = 'hyzd_xz' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hyzd,'xz')
  
    char3 = 'hzzd_xy' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hzzd,'xy')
    char3 = 'hzzd_xz' // trim(char2) // '.txt'
    call IO_output_2D_general(char3, 'g', 'g', hzzd,'xz')
  end if 
!
  deallocate(aps)
!
end subroutine IO_output_2D_xyz
