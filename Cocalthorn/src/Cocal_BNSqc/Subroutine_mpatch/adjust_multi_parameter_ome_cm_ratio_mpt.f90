subroutine adjust_multi_parameter_ome_cm_ratio_mpt(flag_param,istep_niq)
  use phys_constant, only  : long
  use grid_parameter, only : mass_eps, eps
  use interface_minv
  use interface_adjust_copy_ome_cm_ratio_from_mpt
  use interface_adjust_copy_ome_cm_ratio_to_mpt
  use interface_adjust_calc_fncval_Virial_Py_Mratio_mpt
  use interface_error_adjust_parameter
  implicit none
  integer, parameter :: niq = 3
  integer :: flag_param, istep_niq
  integer :: ii, i, j, k
  integer, save :: count_adj
  real(long) :: fac0, error_param
  real(long), save :: msec_x_oold(niq), msec_x_old(niq)
  real(long), save :: msec_f_oold(niq), msec_f_old(niq)
  real(long), save :: msec_dx(niq), msec_x_der(niq), msec_f_der(niq)
  real(long), save :: msec_ff(niq,niq)
  real(long) :: jacobian_at_x_old(niq+1,niq+1), rhs_at_x_old(niq)
  real(long) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
!
  if (istep_niq.eq.-1) then 
    call adjust_copy_ome_cm_ratio_from_mpt(niq,msec_x_oold)
    call adjust_calc_fncval_Virial_Py_Mratio_mpt(niq,msec_f_oold)
    msec_x_old(1:niq) = msec_x_oold(1:niq) + 0.05*msec_x_oold(1:niq)
    call adjust_copy_ome_cm_ratio_to_mpt(niq,msec_x_old)
    istep_niq = istep_niq + 1
    flag_param = 1
    count_adj = 0
    return
  end if
!
!write(6,*)'istep_niq', istep_niq
!write(6,'(a5,1p,3e14.6)')'xoold', (msec_x_oold(ii),ii=1,niq)
!write(6,'(a5,1p,3e14.6)')'x_old', (msec_x_old(ii),ii=1,niq)
!write(6,'(a5,1p,3e14.6)')'f_old', (msec_f_old(ii),ii=1,niq)
!write(6,'(a5,1p,3e14.6)')'x_der', (msec_x_der(ii),ii=1,niq)
!write(6,'(a5,1p,3e14.6)')'f_der', (msec_f_der(ii),ii=1,niq)
!!
  if (istep_niq.eq.0) then 
    call adjust_calc_fncval_Virial_Py_Mratio_mpt(niq,msec_f_old)
  end if
  if (istep_niq.ge.1) then
    ii = istep_niq
    call adjust_calc_fncval_Virial_Py_Mratio_mpt(niq,msec_f_der)
    msec_dx(ii)       = msec_x_oold(ii) - msec_x_old(ii)     ! x_{n-1}-x_n
    msec_ff(1:niq,ii) = msec_f_der(1:niq)
  end if
!
! Set next parameters
  if (istep_niq.ge.0.and.istep_niq.lt.niq) then
    ii = istep_niq + 1
    msec_x_der(1:niq) = msec_x_old(1:niq)
    msec_x_der(ii)    = msec_x_oold(ii)      ! off diagonal point
    call adjust_copy_ome_cm_ratio_to_mpt(niq,msec_x_der)
    istep_niq = istep_niq + 1
    flag_param = 1
    return
  end if
!
! computation of Jacobian at point x1 i.e Xold
  do i = 1, niq
    do j = 1, niq
      jacobian_at_x_old(i,j) = (msec_ff(i,j) - msec_f_old(i))/msec_dx(j)
    end do
  end do
!   computation of the RHS: A*x_n - F_n
  do j = 1, niq
    rhs_at_x_old(j) = 0.0d0
    do k = 1, niq
      rhs_at_x_old(j) = rhs_at_x_old(j) + jacobian_at_x_old(j,k)*msec_x_old(k)
    end do
    rhs_at_x_old(j) = rhs_at_x_old(j) - msec_f_old(j)
  end do
!
  call minv(jacobian_at_x_old, rhs_at_x_old, niq, niq+1)  ! at the end B=Xnew
!
! update parameters
  count_adj = count_adj + 1
  ii = min0(5,count_adj)
  fac0 = facp(ii)
  msec_x_oold(1:niq) = msec_x_old(1:niq)
  msec_x_old( 1:niq) = fac0*rhs_at_x_old(1:niq)+(1.0d0-fac0)*msec_x_oold(1:niq)
  call adjust_copy_ome_cm_ratio_to_mpt(niq,msec_x_old)
  istep_niq = 0
!
  call error_adjust_parameter(niq,msec_x_old,msec_x_oold, &
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
end subroutine adjust_multi_parameter_ome_cm_ratio_mpt
