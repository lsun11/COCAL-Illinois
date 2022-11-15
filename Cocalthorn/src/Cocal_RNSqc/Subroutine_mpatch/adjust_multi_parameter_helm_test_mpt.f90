subroutine adjust_multi_parameter_helm_test_mpt(flag_param,istep_niq)
  use phys_constant, only  : long
  use grid_parameter, only : mass_eps, eps
  use interface_minv
  use interface_adjust_copy_helm_test_from_mpt
  use interface_adjust_copy_helm_test_to_mpt
  use interface_adjust_calc_fncval_helm_test_mpt
  use interface_error_adjust_parameter
  implicit none
  integer, parameter :: niq = 1
  integer :: flag_param, istep_niq
  integer :: ii, i, j, k
  integer, save :: count_adj
  real(long) :: fac0, error_param
  real(long), save :: msec_x_oold(niq), msec_x_old(niq)
  real(long), save :: msec_f_oold(niq), msec_f_old(niq)
  real(long), save :: msec_dx(niq), msec_x_der(niq), msec_f_der(niq)
  real(long), save :: msec_ff(niq,niq)
  real(long) :: jacobian_at_x_old(niq+1,niq+1), rhs_at_x_old(niq)
  real(long) :: jacobian
!m  real(long) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
  real(long) :: facp(5) = (/ 2.0d-1, 5.0d-1, 1.0d-0, 1.0d-0, 1.0d-0 /)
!
  if (istep_niq.eq.-1) then 
    call adjust_copy_helm_test_from_mpt(niq,msec_x_oold)
    call adjust_calc_fncval_helm_test_mpt(niq,msec_f_oold)
    msec_x_old(1:niq) = msec_x_oold(1:niq) + 0.05*msec_x_oold(1:niq)
    call adjust_copy_helm_test_to_mpt(niq,msec_x_old)
    istep_niq = istep_niq + 1
    flag_param = 1
    count_adj = 0
    return
  end if
!
  i = 1
  call adjust_calc_fncval_helm_test_mpt(niq,msec_f_old)
  msec_dx(i) = msec_x_old(i) - msec_x_oold(i)     ! x_{n-1}-x_{n-2}
  jacobian   =(msec_f_old(i) - msec_f_oold(i))/msec_dx(i)
  rhs_at_x_old(1:niq) = msec_x_old(1:niq) - msec_f_old(1:niq)/jacobian
!
    write(6,'(1a20,1p,2e12.4)')' -- OLD -> NEW  --  ', &
    &                          msec_x_oold(i), msec_x_old(i)
    write(6,'(1a20,1p,2e12.4)')' -- OLD -> NEW  --  ', &
    &                          msec_f_oold(i), msec_f_old(i)
!   computation of the RHS: A*x_n - F_n
! update parameters
  count_adj = count_adj + 1
  ii = min0(5,count_adj)
  fac0 = facp(ii)
  msec_x_der(1:niq)  = msec_x_oold(1:niq)
  msec_x_oold(1:niq) = msec_x_old(1:niq)
  msec_f_oold(1:niq) = msec_f_old(1:niq)
  msec_x_old( 1:niq) = fac0*rhs_at_x_old(1:niq)+(1.0d0-fac0)*msec_x_old(1:niq)
  call adjust_copy_helm_test_to_mpt(niq,msec_x_old)
  istep_niq = 0
!
  call error_adjust_parameter(niq,msec_x_old,msec_x_der, &
  &                           error_param,flag_param)
!
  write(6,'(1a30,2(7x,i5))') ' -- parameter iteration # --  ',count_adj
  do ii = 1, niq
    write(6,'(1a20,1p,2e12.4)')' -- OLD -> NEW  --  ', &
    &                          msec_x_oold(ii), msec_x_old(ii)
  end do
!  write(6,'(1a20,1p,2e12.4)')' -- ADM and Komar-- ', &
!  &                          admmass_asymp,  komarmass_asymp
!  write(6,'(1a20,1p,2e12.4)')' -- 1 - MK/MADM  -- ', fnc_val_bak, fnc_val
!
end subroutine adjust_multi_parameter_helm_test_mpt
