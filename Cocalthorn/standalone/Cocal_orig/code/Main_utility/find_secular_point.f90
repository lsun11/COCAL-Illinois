!______________________________________________
include '../Module/phys_constant.f90'
include '../Module/def_quantities.f90'
include '../Module/def_quantities_derived.f90'
include '../Subroutine/IO_input_physq_plot.f90'
include '../Subroutine/IO_output_physq_plot.f90'
include '../Analysis/Subroutine/IO_output_physq_paper.f90'
include '../Analysis/Module/def_physq_to_array.f90'
!
!______________________________________________
!
!       pick up data in a sequence
!______________________________________________
PROGRAM find_secular_point
!
  use def_physq_to_array
  implicit none
  integer :: flag_restmass = 0
  integer :: total_num_sol, nn = 100
  integer :: ii, jj, num_physq, iseq_end, iseq_sig
  real(8) :: katamuki, param
!
  do ii = 1, nn
    jj = ii
    call IO_input_physq_plot(jj,flag_restmass)
    if (flag_restmass.eq.9999) then 
      total_num_sol = ii - 1 
      exit
    end if
    call physq_to_array(ii,flag_restmass,num_physq)
  end do
!
  if (physq_array(9,1).gt.physq_array(9,2)) then
    write(6,*) ' Axis ratio decreasing order ' 
    iseq_end = 0
    iseq_sig = 1
  else 
    write(6,*) ' Axis ratio increasing order ' 
    iseq_end = total_num_sol + 1
    iseq_sig = - 1
  end if
!
  physq_array(:,iseq_end) = 0.0d0
  read(5,*) param
!
  do jj = 3, num_physq
    katamuki = (physq_array(jj,iseq_end+iseq_sig*2)  & 
  &           - physq_array(jj,iseq_end+iseq_sig*3)) & 
  &           /(physq_array( 9,iseq_end+iseq_sig*2)  &
  &           - physq_array( 9,iseq_end+iseq_sig*3))
    physq_array(jj,iseq_end) = katamuki              &
!  &   *(1.0d0 - physq_array( 9,iseq_end+iseq_sig*1)) & 
  &   *(param - physq_array( 9,iseq_end+iseq_sig*1)) & 
  &           + physq_array(jj,iseq_end+iseq_sig*1)
  end do
!
  call array_to_physq(iseq_end+iseq_sig*1,flag_restmass,num_physq)
  call IO_output_physq_plot(1,flag_restmass)
  call array_to_physq(iseq_end,flag_restmass,num_physq)
  call IO_output_physq_plot(2,flag_restmass)
!
  call IO_output_physq_paper(1)
!
END PROGRAM find_secular_point
