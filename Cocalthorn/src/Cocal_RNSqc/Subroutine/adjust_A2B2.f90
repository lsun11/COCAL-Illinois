subroutine adjust_A2B2(flag_st,istep_st)
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntfeq, npfxzp
  use def_matter, only : omef
  use def_matter_parameter, only : ome, A2DR, DRAT_A2DR, B2DR, DRAT_B2DR
  use interface_minv
  use make_array_1d
  use interface_search_max_radial
  implicit none
  integer, parameter :: nst = 2
  integer :: flag_st, istep_st
  integer :: ii, i, j, k
  integer, save :: count_st
  real(long) :: fac0, error_st
  real(long), save :: msec_x_oold(nst), msec_x_old(nst)
  real(long), save :: msec_f_oold(nst), msec_f_old(nst)
  real(long), save :: msec_dx(nst), msec_x_der(nst), msec_f_der(nst)
  real(long), save :: msec_ff(nst,nst)
  real(long) :: jacobian_at_x_old(nst+1,nst+1), rhs_at_x_old(nst)
  real(long) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
!facp  real(long) :: facp(8) = (/         2.0d-2, 4.0d-2, 7.0d-2,        &
!facp  &                          1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
  real(long), pointer :: omef_x(:)
  real(long)          :: omef_max, x_omax, omef_eq
  integer             :: irf_omax
!
! -- secant method: iteration for A2DR and B2DR
!
  call alloc_array1d(omef_x,0,nrf)
  omef_x(0:nrf)=omef(0:nrf,ntfeq,npfxzp)
  irf_omax = maxloc(omef_x,DIM=1) - 1
!  write(6,*) maxloc(omef_x) - 1
!  write(6,*)'ome max',omef_x(irf_omax-1),omef_x(irf_omax),omef_x(irf_omax+1)
!stop
  call search_max_radial(irf_omax,omef_x,omef_max,x_omax)
  omef_eq = omef(nrf,ntfeq,npfxzp)
  deallocate(omef_x)
!
  if (istep_st.eq.-1) then
    msec_x_oold(1) = A2DR
    msec_x_oold(2) = B2DR
    msec_f_oold(1) = omef_eq /ome - DRAT_A2DR
    msec_f_oold(2) = omef_max/ome - DRAT_B2DR
    msec_x_old(1:nst) = msec_x_oold(1:nst) + 0.01*msec_x_oold(1:nst)
    A2DR = msec_x_old(1)
    B2DR = msec_x_old(2)
    istep_st = istep_st + 1
    flag_st  = 1
    count_st = 0
    return
  end if
!
write(6,*)'istep A2B2', istep_st
write(6,'(a5,1p,3e14.6)')'xoold', (msec_x_oold(ii),ii=1,nst)
write(6,'(a5,1p,3e14.6)')'x_old', (msec_x_old(ii),ii=1,nst)
write(6,'(a5,1p,3e14.6)')'f_old', (msec_f_old(ii),ii=1,nst)
write(6,'(a5,1p,3e14.6)')'x_der', (msec_x_der(ii),ii=1,nst)
write(6,'(a5,1p,3e14.6)')'f_der', (msec_f_der(ii),ii=1,nst)
!
  if (istep_st.eq.0) then
    msec_f_old(1) = omef_eq /ome - DRAT_A2DR
    msec_f_old(2) = omef_max/ome - DRAT_B2DR
  end if
  if (istep_st.ge.1) then
    ii = istep_st
    msec_f_der(1) = omef_eq /ome - DRAT_A2DR
    msec_f_der(2) = omef_max/ome - DRAT_B2DR
    msec_dx(ii)       = msec_x_oold(ii) - msec_x_old(ii)  ! x_{n-1}-x_n     
    msec_ff(1:nst,ii) = msec_f_der(1:nst)  
  end if
!
! Set next parameters
  if (istep_st.ge.0.and.istep_st.lt.nst) then
    ii = istep_st + 1
    msec_x_der(1:nst) = msec_x_old(1:nst)
    msec_x_der(ii)    = msec_x_oold(ii)      ! off diagonal point
    A2DR = msec_x_der(1)
    B2DR = msec_x_der(2)
    istep_st = istep_st + 1
    flag_st = 1
    return
  end if
!
!
! computation of Jacobian at point x1 i.e Xold
  do i = 1, nst
    do j = 1, nst
      jacobian_at_x_old(i,j) = (msec_ff(i,j) - msec_f_old(i))/msec_dx(j)
    end do
  end do
!   computation of the RHS: A*x_n - F_n
  do j = 1, nst
    rhs_at_x_old(j) = 0.0d0
    do k = 1, nst
      rhs_at_x_old(j) = rhs_at_x_old(j) + jacobian_at_x_old(j,k)*msec_x_old(k)
    end do
    rhs_at_x_old(j) = rhs_at_x_old(j) - msec_f_old(j)
  end do
!
  call minv(jacobian_at_x_old, rhs_at_x_old, nst, nst+1)  ! at the end B=Xnew
!
! update parameters
  count_st = count_st + 1
  ii = min0(5,count_st)
!facp  ii = min0(8,count_st)
  fac0 = facp(ii)
  msec_x_oold(1:nst) = msec_x_old(1:nst)
  msec_x_old( 1:nst) = fac0*rhs_at_x_old(1:nst)+(1.0d0-fac0) &
  &                        *msec_x_oold(1:nst)
  A2DR = msec_x_old(1)
  B2DR = msec_x_old(2)
  istep_st = 0
!                                                                             
  call error_adjust_parameter(nst,msec_x_old,msec_x_oold, &
  &                           error_st,flag_st)
!
  write(6,'(1a30,2(7x,i5))') ' -- parameter iteration # --  ',count_st
  write(6,'(1a30,1p,2e12.4)')' --   A2DR    OLD -> NEW  --  ', &
  &                            msec_x_oold(1), msec_x_old(1)
  write(6,'(1a30,1p,2e12.4)')' --   B2DR    OLD -> NEW  --  ', &
  &                            msec_x_oold(2), msec_x_old(2)
!
  return
end subroutine adjust_A2B2
