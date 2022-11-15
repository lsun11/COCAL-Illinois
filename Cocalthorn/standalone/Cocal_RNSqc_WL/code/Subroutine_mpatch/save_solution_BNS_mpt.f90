subroutine save_solution_BNS_mpt(icycle)
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

  call write_last_physq_BNS_mpt
  char3 = 'main_bnsphys_all_mpt.txt'
  call printout_physq_BNS_all_mpt(char3)

  call chdir('../')

end subroutine save_solution_BNS_mpt
