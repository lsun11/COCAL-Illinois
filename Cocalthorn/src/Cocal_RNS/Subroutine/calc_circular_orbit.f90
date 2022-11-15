subroutine calc_circular_orbit(total_iteration)
  use phys_constant, only  : long
  use def_bh_parameter, only : ome_bh
  use def_quantities
  implicit none
  real(long) :: xa,xb,fa,fb,eps,delta,f1,f2
  integer :: total_iteration, iter_count, i, icycle, Nmax=50
  eps = 1.0d-07
  total_iteration = 0
! 
  xa = ome_bh
  call iteration_BBH_CF(iter_count)
  total_iteration = total_iteration + iter_count
  write(6,'(a21,i5)') '--- total iteration =', total_iteration
  call calc_physical_quantities_BBH_CF
  call save_solution(0)
  call write_last_physq
  call printout_physq_console_BBH
  fa = (admmass - komarmass)/komarmass
  ome_bh =ome_bh + 0.2d0*ome_bh 
!
  xb = ome_bh
  call iteration_BBH_CF(iter_count)
  total_iteration = total_iteration + iter_count
  write(6,'(a21,i5)') '--- total iteration =', total_iteration
  call calc_physical_quantities_BBH_CF
  call save_solution(1)
  call write_last_physq
  call printout_physq_console_BBH
  fb = (admmass - komarmass)/komarmass
!
  write(6,*) '--------- Initial values for Omega and f(Omega)=(admmass-komarmass)/komarmass ---------------------------'
  write(6,'(a7,i2,1p,2e20.12)') 'cycle =',0,xa,fa
  write(6,'(a7,i2,1p,2e20.12)') 'cycle =',1,xb,fb
!
  if ( dabs(fa)>dabs(fb) ) then  ! interchange xa,xb  and fa,fb
    f1 = xa
    xa = xb
    xb = f1

    f1 = fa
    fa = fb
    fb = f1
  endif
!
  icycle = 1
  do i = 2,Nmax
    icycle = icycle + 1
    if ( dabs(fa)>dabs(fb) ) then
      f1 = xa
      xa = xb
      xb = f1
      f1 = fa
      fa = fb
      fb = f1
    end if
    delta = (xb-xa)/(fb-fa)
    xb = xa
    fb = fa
    delta = delta*fa
    if ( dabs(delta) < eps ) then
      write(6,*) 'Convergence'
      write(6,*) '------------------------------------------------------------------'
      write(6,*) 'Omega orbit = ',ome_bh
      call printout_physq_console_BBH
      return
    end if
    if (i.lt.10) then
      xa = xa - i*delta/10.0d0
    else
      xa = xa - delta
    end if 
    ome_bh = xa
    write(6,*) '------------------------------------------------------------------'
    write (6,'(a7,i2,a20,e20.12)') 'cycle =', icycle, '        New Omega = ',ome_bh
    call iteration_BBH_CF(iter_count)
    total_iteration = total_iteration + iter_count
    write(6,'(a21,i5)') '--- total iteration =', total_iteration
    call calc_physical_quantities_BBH_CF
    call save_solution(icycle)
    call write_last_physq
    fa = (admmass - komarmass)/komarmass
    write(6,'(a7,i2,1p,2e20.12)')  'cycle =',icycle,xa,fa
    call printout_physq_console_BBH
  end do
  write(6,*) 'No convergence after ',Nmax,' iterations.'

end subroutine calc_circular_orbit
