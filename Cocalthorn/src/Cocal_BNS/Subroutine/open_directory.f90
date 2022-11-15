subroutine open_directory(ia,ib)
  use phys_constant, only  : long
  use def_quantities
  implicit none 
  integer :: ia, ib
  character(40) :: char1, char2, char3, char4, char5
  character(100) :: dircommand

  write(char1, '(i5)') ia
  char2 = adjustl(char1)
  if (ia < 10) then
    char5 = 'sol_0' // trim(char2) 
  else
    char5 = 'sol_' // trim(char2)
  endif

  dircommand = 'mkdir ' // char5
  call system(dircommand)

!  The following doesn't work *********************************
!  dircommand = 'cd ' // char5
!  call system(dircommand)
  call chdir(char5)

!  call system('pwd')
!  call chdir('../')

end subroutine open_directory
