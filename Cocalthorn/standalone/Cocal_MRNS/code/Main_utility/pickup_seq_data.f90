!______________________________________________
include '../Module/phys_constant.f90'
include '../Module/def_quantities.f90'
include '../Module/def_quantities_derived.f90'
include '../Subroutine/IO_input_physq_plot.f90'
include '../Subroutine/IO_output_physq_plot.f90'
include '../Analysis/Subroutine/IO_output_physq_paper.f90'
!
!______________________________________________
!
!       pick up data in a sequence
!______________________________________________
PROGRAM pick_up_seq_data
!
  implicit none
  integer :: flag_restmass
  integer :: total_num_sol, pickup_num_sol, sol_id(40)
  integer :: ii, jj
!
  open(20,file='rns_pickup_seq.dat',status='old')
  read(20,'(2i5)') total_num_sol, pickup_num_sol
  do ii = 1, pickup_num_sol
    read(20,'(i5)') sol_id(ii)
  end do
!
  do jj = 1, total_num_sol
    call IO_input_physq_plot(jj,flag_restmass)
    do ii = 1, pickup_num_sol
      if (jj.eq.sol_id(ii)) then 
        call IO_output_physq_plot(ii,flag_restmass)
        call IO_output_physq_paper(ii)
      end if
    end do
  end do
!
END PROGRAM pick_up_seq_data
