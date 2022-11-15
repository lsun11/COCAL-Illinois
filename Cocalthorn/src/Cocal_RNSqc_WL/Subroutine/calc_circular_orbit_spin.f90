subroutine calc_circular_orbit_spin(total_iteration,niq)
  use phys_constant, only  : long
  use def_bh_parameter, only : ome_bh, spin_bh
  use def_quantities
  use def_iter_quantities
  use make_array_1d
  use make_array_2d
  use interface_msec_store_f_vector
  use interface_msec_copy_to_iter_quantities
  use interface_minv
  implicit none
  real(long) :: delta, cf, edis

  real(long), pointer :: jacobian_at_x_old(:,:), msec_ff(:,:)
  real(long), pointer :: rhs_at_x_old(:), msec_dx(:), msec_x_der(:), msec_f_der(:)

  integer :: total_iteration, iter_count, Nmax=30
  integer :: istep, ii, i, j, k, niq

  call alloc_array2d(jacobian_at_x_old,1,niq+1,1,niq+1)
  call alloc_array2d(msec_ff          ,1,niq  ,1,niq)
  call alloc_array1d(rhs_at_x_old,1,niq)
  call alloc_array1d(msec_dx     ,1,niq)
  call alloc_array1d(msec_x_der  ,1,niq)
  call alloc_array1d(msec_f_der  ,1,niq)

  delta = 1.0d-10
!
  istep = 0
  write(6,'(a13,i5,1p,10e20.12)')  'First point: ', istep, (msec_x_oold(i),i=1,niq)
  if (niq.eq.1) then
    call msec_copy_to_iter_quantities(msec_x_oold)  ! It is initialized in allocate_iter_quantities
    call iteration_BBH_CF(iter_count)
    total_iteration = total_iteration + iter_count
    write(6,'(a21,i5)') '--- total iteration =', total_iteration
    call calc_physical_quantities_BBH_CF
    call msec_store_f_vector(msec_f_oold)  
    call save_solution(0)
    call write_last_physq
    call printout_physq_console_BBH
  endif
  write(6,*) '-----------------------------------------------------------------------------------------------'
!
  do i=1,niq
    msec_x_old(i) = msec_x_oold(i) + 0.05*msec_x_oold(i)
  end do
!
  istep = 1
  write(6,'(a11,i5,1p,10e20.12)')  'New point: ', istep, (msec_x_old(i),i=1,niq)
  call msec_copy_to_iter_quantities(msec_x_old)
  call iteration_BBH_CF(iter_count)
  total_iteration = total_iteration + iter_count
  write(6,'(a21,i5)') '--- total iteration =', total_iteration
  call calc_physical_quantities_BBH_CF
  call msec_store_f_vector(msec_f_old)  
  call save_solution(1)
  call write_last_physq
  call printout_physq_console_BBH
  write(6,*) '-----------------------------------------------------------------------------------------------'
  
  do ii = 2,Nmax
    istep = istep+1
    if (niq > 1) then
      do i=1,niq
        msec_dx(i) = msec_x_oold(i) - msec_x_old(i)  ! x_{n-1}-x_n
        msec_x_der(1:niq) = msec_x_old(1:niq)
        msec_x_der(i) = msec_x_oold(i)        ! off diagonal point used for the calc of derivative
!
        write(6,'(a11,i5,i5,1p,10e20.12)') 'Add point: ', istep, i, (msec_x_der(j),j=1,niq)
        call msec_copy_to_iter_quantities(msec_x_der)   ! store new values of omega,spin,...
        call iteration_BBH_CF(iter_count)
        total_iteration = total_iteration + iter_count
        write(6,'(a21,i5)') '--- total iteration =', total_iteration
        call calc_physical_quantities_BBH_CF
        call msec_store_f_vector(msec_f_der)  ! store values of quantities to be solved with secant method
!        call save_solution(1)
        call write_last_physq
        call printout_physq_console_BBH
        write(6,*) '************************************************************************************************'
!
        msec_ff(1:niq,i) = msec_f_der(1:niq)
      end do
    else
      msec_dx(1)   = msec_x_oold(1) - msec_x_old(1)
      msec_ff(1,1) = msec_f_oold(1)
    endif
!   computation of Jacobian at point x1 i.e Xold
    do i=1,niq
      do j=1,niq
        jacobian_at_x_old(i,j) = (msec_ff(i,j) - msec_f_old(i))/msec_dx(j)
      end do
    end do
!   computation of the RHS: A*x_n - F_n
    do j=1,niq
      rhs_at_x_old(j) = 0.0d0
      do k=1,niq
        rhs_at_x_old(j) = rhs_at_x_old(j) + jacobian_at_x_old(j,k)*msec_x_old(k)
      end do
      rhs_at_x_old(j) = rhs_at_x_old(j) - msec_f_old(j)
    end do
!
    call minv(jacobian_at_x_old, rhs_at_x_old, niq, niq+1)  ! at the end B=Xnew
!
    msec_x_oold(1:niq) = msec_x_old(1:niq)
    msec_f_oold(1:niq) = msec_f_old(1:niq)   ! this is not necessary actually. msec_f_oold is not needed any more 

!   the newly calculated point
    if (ii.lt.10) then
      cf = dble(ii)/10.0d0
      msec_x_old(1:niq) = cf*rhs_at_x_old(1:niq) + (1.0d0-cf)*msec_x_oold(1:niq)
    else
      msec_x_old(1:niq) = rhs_at_x_old(1:niq) 
    end if

    write(6,'(a11,i5,1p,10e20.12)')  'New point: ', istep, (msec_x_old(i),i=1,niq)
    call msec_copy_to_iter_quantities(msec_x_old)
    call iteration_BBH_CF(iter_count)
    total_iteration = total_iteration + iter_count
    write(6,'(a21,i5)') '--- total iteration =', total_iteration
    call calc_physical_quantities_BBH_CF
    call msec_store_f_vector(msec_f_old)
    call save_solution(istep)
    call write_last_physq
    call printout_physq_console_BBH
    write(6,*) '-----------------------------------------------------------------------------------------------'

!   check the distance between msec_x_oold and msec_x_old
    edis = 0
    do i=1,niq
      edis = edis + (msec_x_oold(i) - msec_x_old(i))*(msec_x_oold(i) - msec_x_old(i))
    end do
    edis = sqrt(edis)
    if (edis < delta) then
      write(6,*) '================================================================================================'
      write(6,'(a13,i5,1p,10e20.12)') 'Convergence: ', istep, (msec_x_old(i),i=1,niq)
      exit
    end if
  end do
  if (istep.eq.Nmax) then
    if (edis.lt.delta)  then
      write(6,*) 'Convergence after ',istep-1,' iterations.'
    else
      write(6,*) 'No convergence after ',Nmax-1,' iterations.'
    endif
  else
     write(6,*) 'Convergence after ',istep-1,' iterations.'
  endif
  write(6,*) 'Deallocating jacobian_at_x_old,ff,rhs_at_x_old...'
  deallocate(jacobian_at_x_old, msec_ff, rhs_at_x_old, msec_dx, msec_x_der, msec_f_der)
end subroutine calc_circular_orbit_spin
