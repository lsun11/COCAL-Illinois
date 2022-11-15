subroutine save_solution(icycle)
  use phys_constant, only  : long
  use def_bh_parameter
  use def_quantities
  use grid_parameter, only : rgin
  implicit none 
  integer :: icycle
  character(40) :: char1, char2, char3, char4, char5
  character(100) :: dircommand

  write(char1, '(i5)') icycle
  char2 = adjustl(char1)
!  write(char3, '(e20.12)') omega_bh
!  char4 = adjustl(char3)
!  char5 = 'cycle_' // trim(char2) // '_omega_'// trim(char4)
  if (icycle < 10) then
    char5 = 'cycle_0' // trim(char2) 
  else
    char5 = 'cycle_' // trim(char2)
  endif

  dircommand = 'mkdir ' // char5
  call system(dircommand)

!  The following doesn't work *********************************
!  dircommand = 'cd ' // char5
!  call system(dircommand)

  call chdir(char5)
!  call system('pwd')
!  call IO_output_solution_3D
!  call IO_output_BBH_CF    ! for plot_x.dat,...

  open(15,file='global_' // trim(char5) // '.txt',status='unknown')
  write(15,'(1p,13e16.8)') rgin, spin_bh, ome_bh, admmass, komarmass,admmass_thr,  &
                   &       app_hor_area_bh, irredmass, bindingene,   &
                   &       angmom, angmom_thr, angmom_smarr, qua_loc_spin
  close(15)
  call chdir('../')

end subroutine save_solution
