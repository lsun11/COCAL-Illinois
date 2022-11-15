subroutine next_solution_BNS_mpt(iseq)
  use phys_constant,  only : long
  use def_peos_parameter
  use def_peos_parameter_mpt
  use def_matter_parameter_mpt
  use def_matter_parameter
  implicit none
  integer :: impt, iseq
  character(len=1) :: np(2) = (/'1', '2'/)
  character(LEN=200) :: inputfile, cmd, tempchar1, tempchar2
!
  do impt = 1, 2
    call copy_def_matter_parameter_from_mpt(impt)
    emdc = emdc + 0.05d0*emdc
    call copy_def_matter_parameter_to_mpt(impt)
!
!    Changing the initial central density for TOV 
!    call copy_def_peos_parameter_from_mpt(impt)
!    rhoini_cgs = rhoini_cgs_iseq + 1.0d+16*dble(iseq)   
!    write(tempchar1, '(es13.5)')  rhoini_cgs
!    cmd = " sed -i ""1s:\(.\{13\}\)\(.\{13\}\):\1" // trim(tempchar1) // ":"" peos_parameter_mpt" //np(impt)// ".dat"
!    call system(cmd)
!    call peos_initialize_mpt(impt)
!    call copy_def_peos_parameter_to_mpt(impt)
  end do
!
end subroutine next_solution_BNS_mpt
