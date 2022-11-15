!______________________________________________
include '../Module/phys_constant.f90'
include '../Module/def_quantities.f90'
include '../Module/def_quantities_derived.f90'
include '../Subroutine/IO_input_physq_plot.f90'
include '../Analysis/Module/def_physq_to_array.f90'
!
!______________________________________________
!
!       pick up data in a sequence
!______________________________________________
PROGRAM take_data
!
  use def_physq_to_array
  implicit none
  integer :: flag_restmass = 0
  integer :: total_num_sol, nn = 100
  integer :: ii, jj, num_physq
  integer :: param(100), num_output_data
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
  read(5,*) num_output_data, (param(ii), ii = 1, num_output_data)
!
  open(1,file='data.dat',status='old')
  write(1,'(1p,9es18.10)')(physq_array(param(ii),1), ii = 1, num_output_data)
  close(1)
!
END PROGRAM take_data
