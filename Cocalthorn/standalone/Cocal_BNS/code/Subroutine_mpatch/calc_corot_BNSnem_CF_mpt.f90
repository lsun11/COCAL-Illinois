subroutine calc_corot_BNSnem_CF_mpt(total_iteration, sw_master)
  use phys_constant, only  : long, nmpt
  use def_iter_quantities
  use grid_parameter
  use def_quantities
  use def_matter_parameter
  use def_matter_parameter_mpt
  use interface_copy_msec_BNS_iterqt_to_mpt
  use interface_copy_msec_BNS_iterqt_from_mpt
  use interface_store_msec_BNS_fiterqt_mpt
  use interface_violation_gridpoint_MoC_CF_peos_corot
  use interface_violation_gridpoint_HaC_CF_peos_corot
  use interface_IO_output_3D_general
  use interface_IO_output_2D_general
  use interface_IO_output_1D_general
  use make_array_1d
  use make_array_2d
  use make_array_3d
  use make_array_4d
  use interface_minv
  implicit none
  real(long), pointer :: pot(:,:,:), HaC_vio(:,:,:), MoC_vio(:,:,:,:)
  real(long) :: xa,xb,fa,fb,rm_eps,delta,f1,f2, iter_eps
  real(long) :: cf, edis_iterqt, maxfold
  real(long) :: cycleeps(7) = (/ 1.0d0, 1.0d-1, 1.0d-1, 1.0d-2, 1.0d-3, 1.0d-4, 1.0d-5 /)
  integer :: total_iteration, iter_count, i, ii, icycle, j, k, Nmax=20
  integer :: impt, niq
  character(3) :: sw_master
  character(30) :: char1, char2, char3, char4, char5
  real(long), pointer :: jacobian_at_x_old(:,:), msec_ff(:,:)
  real(long), pointer :: rhs_at_x_old(:), msec_dx(:), msec_x_der(:), msec_f_der(:)

  select case (sw_master)
    case ("010", "030", "050")  ! iteration over one restmass or adm mass, or compactness 
      niq = 1                   ! emdc1=emdc2
    case ("020", "040", "060")  ! iteration over two restmasses or adm masses, or compactnesses 
      niq = 4                   ! emdc1, emdc2, cm, rs2
    case ("110", "130", "150", "210", "230", "250") ! iteration over distance and one mass
      niq = 2                                       ! rs1=rs2, emdc1=emdc2
    case ("120", "140", "160", "220", "240", "260") ! iteration over distance and two masses
      niq = 5                                       ! rs1, emdc1, emdc2, cm, rs2
    case default
      write(6,*) "Check switches sw_sepa, sw_quant, sw_spin. Unrecognized options",sw_master,"...exiting."
      stop
  end select

  write(6,*) "Master switch, niq =====>", sw_master, " -- ", niq

  call allocate_BNS_iter_quantities(niq)

  call alloc_array2d(jacobian_at_x_old,1,niq+1,1,niq+1)
  call alloc_array2d(msec_ff          ,1,niq  ,1,niq)
  call alloc_array1d(rhs_at_x_old,1,niq)
  call alloc_array1d(msec_dx     ,1,niq)
  call alloc_array1d(msec_x_der  ,1,niq)
  call alloc_array1d(msec_f_der  ,1,niq)

  delta = 1.0d-06
  rm_eps = 1.0d-06
  total_iteration = 0
  char3 = 'bnsphyseq.dat'
!
  icycle = 0
  write(6,'(a53,i4,a53)')   "------------------------------------------- Cycle", icycle, &
                       &    "-------------------------------------------------"
  call copy_msec_BNS_iterqt_from_mpt(msec_x_oold, sw_master)  ! Initializes msec_x_oold
  write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  First point    : ', (msec_x_oold(i),i=1,niq)
  if (niq.eq.1) then
    call open_directory_mpt(0)
    iter_eps = cycleeps(1)
    call iter_corot_BNS_CF_mpt(iter_count,1,iter_eps)         !  1 for freezing hydro in first 5 iterations
    total_iteration = total_iteration + iter_count
    write(6,'(a11,i2,a25,i4,a32,i5)') '--- cycle =', icycle, '   number of iterations =', iter_count,  &
        &  '   total number of iterations = ', total_iteration
    call calc_physical_quantities_BNS_CF_mpt   
    call store_msec_BNS_fiterqt_mpt(msec_f_oold, sw_master)         ! Initializes msec_f_oold
    write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  f(First point) : ', (msec_f_oold(i),i=1,niq)
    call write_last_physq_BNS_mpt
    call printout_physq_BNS_all_mpt(char3)
    call chdir('../')
  endif
!
  do i=1,niq
    msec_x_old(i) = msec_x_oold(i) - 0.01*msec_x_oold(i)
  end do
  call copy_msec_BNS_iterqt_to_mpt(msec_x_old, sw_master)
  if (sw_sepa.ne.0 .or. mod(sw_quant,2).eq.0) then     !  iteration over distance/omega or cm
    write(6,*)  "*** sw_quant=", sw_quant  
    call update_coordinates_BNS_mpt
  endif
!
  icycle = 1
  write(6,'(a53,i4,a53)')   "------------------------------------------- Cycle", icycle, &
                       &    "-------------------------------------------------"
  write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  Second point   : ',  (msec_x_old(i),i=1,niq)
  call open_directory_mpt(1)
  iter_eps = cycleeps(icycle)
  call iter_corot_BNS_CF_mpt(iter_count,1,iter_eps)     !  1 for freezing hydro in first 5 iterations
  total_iteration = total_iteration + iter_count
  write(6,'(a11,i2,a25,i4,a32,i5)') '--- cycle =', icycle, '   number of iterations =', iter_count,  &
      &  '   total number of iterations = ', total_iteration
  call calc_physical_quantities_BNS_CF_mpt   
  call store_msec_BNS_fiterqt_mpt(msec_f_old, sw_master)
  write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  f(Second point): ',  (msec_f_old(i),i=1,niq)
  call write_last_physq_BNS_mpt
  call printout_physq_BNS_all_mpt(char3)
  call chdir('../')

  do ii = 2, Nmax
    icycle = icycle+1
    write(6,'(a53,i4,a53)')   "------------------------------------------- Cycle", icycle, &
                         &    "-------------------------------------------------"

    iter_eps = rm_eps
    if (ii.le.7)     iter_eps = cycleeps(ii)
    write(6,*)  "===== ii, iter_eps=", ii, iter_eps
!
    if (niq > 1) then
      do i=1,niq
        msec_dx(i) = msec_x_oold(i) - msec_x_old(i)  ! x_{n-1}-x_n
        msec_x_der(1:niq) = msec_x_old(1:niq)
        msec_x_der(i) = msec_x_oold(i)        ! off diagonal point used for the calc of derivative
        call copy_msec_BNS_iterqt_to_mpt(msec_x_der, sw_master)
        if (sw_sepa.ne.0 .or. mod(sw_quant,2).eq.0) then     !  iteration over distance/omega or cm
          call update_coordinates_BNS_mpt
        endif
!
        write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  Add point      : ',  (msec_x_der(j),j=1,niq)
        call iter_corot_BNS_CF_mpt(iter_count,icycle,iter_eps)        !  1 for freezing hydro in first 5 iterations
        total_iteration = total_iteration + iter_count
        write(6,'(a11,i2,a25,i4,a32,i5,a8,i1)') '--- cycle =', icycle, '   number of iterations =', iter_count,  &
           &  '   total number of iterations = ', total_iteration, '   niq =', niq
        call calc_physical_quantities_BNS_CF_mpt
        call store_msec_BNS_fiterqt_mpt(msec_f_der, sw_master)         ! Initializes msec_f_der
        write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  f(Add point)   : ',  (msec_f_der(j),j=1,niq)
        call write_last_physq_BNS_mpt
        call printout_physq_BNS_all_mpt(char3)
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

    write(6,*) "______________________________________________________________________________"
    write(6,*) "AAAAA jacob", icycle, "row=1", (jacobian_at_x_old(1,j) ,j=1,niq)
    if (niq==2) write(6,*) "AAAAA jacob", icycle, "row=2", (jacobian_at_x_old(2,j) ,j=1,niq)
    if (niq==3) then
      write(6,*) "AAAAA jacob", icycle, "row=2", (jacobian_at_x_old(2,j) ,j=1,niq)
      write(6,*) "AAAAA jacob", icycle, "row=3", (jacobian_at_x_old(3,j) ,j=1,niq)
    endif

!   computation of the RHS: A*x_n - F_n
    do j=1,niq
      rhs_at_x_old(j) = 0.0d0
      do k=1,niq
        rhs_at_x_old(j) = rhs_at_x_old(j) + jacobian_at_x_old(j,k)*msec_x_old(k)
      end do
      rhs_at_x_old(j) = rhs_at_x_old(j) - msec_f_old(j)
    end do
!
    write(6,*) "AAAAA befor", icycle, (rhs_at_x_old(j),j=1,niq)
    call minv(jacobian_at_x_old, rhs_at_x_old, niq, niq+1)  ! at the end B=Xnew
    write(6,*) "AAAAA after", icycle, (rhs_at_x_old(j),j=1,niq)
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
    call copy_msec_BNS_iterqt_to_mpt(msec_x_old, sw_master)
    if (sw_sepa.ne.0 .or. mod(sw_quant,2).eq.0) then     !  iteration over distance/omega or cm
      call update_coordinates_BNS_mpt
    endif
!
!    if (ii.lt.7) then
!      iter_eps = dmin1(10.0d0**(-ii), edis_iterqt)
!    else
!      iter_eps = rm_eps
!    end if
!
    write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  New point      : ', (msec_x_old(i),i=1,niq)
    call open_directory_mpt(icycle)
    call iter_corot_BNS_CF_mpt(iter_count,icycle,iter_eps)    
    total_iteration = total_iteration + iter_count
    write(6,'(a11,i2,a25,i4,a32,i5)') '--- cycle =', icycle, '   number of iterations =', iter_count,  &
        &  '   total number of iterations = ', total_iteration
    call calc_physical_quantities_BNS_CF_mpt
    call store_msec_BNS_fiterqt_mpt(msec_f_old, sw_master)      
    write(6,'(a11,i2,a19,1p,10e20.12)')  '## icycle =', icycle, '  f(New point)   : ', (msec_f_old(i),i=1,niq)
    call write_last_physq_BNS_mpt
    call printout_physq_BNS_all_mpt(char3)
    call chdir('../')
    write(6,*) '************************************************************************************************'
!
!   Max of relative error wrt the targeted values (except Padm which is absolute error) 
    maxfold = -1.0
    do i=1,niq
      if (maxfold < dabs(msec_f_old(i)) )  maxfold = dabs(msec_f_old(i))
    end do
    write(6,*) "**** icycle, maxfold:", icycle, maxfold
! 
!   Distance between msec_x_oold and msec_x_old
    edis_iterqt = 0
    do i=1,niq
      edis_iterqt = edis_iterqt + (msec_x_oold(i) - msec_x_old(i))*(msec_x_oold(i) - msec_x_old(i))
    end do
    edis_iterqt = sqrt(edis_iterqt)
    write(6,*) "**** icycle, edis_iterqt:", icycle, edis_iterqt
!
!   Convergence when ||x_(i+1)-x_i|| < eps   OR   when max{|f(x_i)|} < eps
    if ((edis_iterqt < rm_eps) .or. (maxfold < rm_eps)) then
      write(6,*) '================================================================================================'
      write(6,'(a11,1p,10e20.12)') '## x_oold: ', (msec_x_oold(i),i=1,niq)
      write(6,'(a11,1p,10e20.12)') '## f_oold: ', (msec_f_oold(i),i=1,niq)
      write(6,'(a22,i2)') '## Convergence cycle: ', icycle
      write(6,'(a11,1p,10e20.12)') '## x_old:  ', (msec_x_old(i),i=1,niq)
      write(6,'(a11,1p,10e20.12)') '## f_old:  ', (msec_f_old(i),i=1,niq)
      exit
    end if
  end do
  if (icycle.eq.Nmax) then
    if (edis_iterqt.lt.rm_eps)  then
      write(6,*) 'Convergence after ',icycle-1,' iterations.'
    else
      write(6,*) 'No convergence after ',Nmax-1,' iterations.'
    endif
  else
     write(6,*) 'Convergence after ',icycle-1,' iterations.'
  endif
 
  write(6,*) 'Deallocating jacobian_at_x_old, msec_ff, rhs_at_x_old, msec_dx, msec_x_der, msec_f_der.'
  deallocate(jacobian_at_x_old, msec_ff, rhs_at_x_old, msec_dx, msec_x_der,msec_f_der)

end subroutine calc_corot_BNSnem_CF_mpt
