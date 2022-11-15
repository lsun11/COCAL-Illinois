subroutine calc_lecc_spin_BNS_CF_mpt(total_iteration)
  use phys_constant, only  : long, nmpt
  use grid_parameter
  use def_quantities
  use def_matter_parameter
  use def_matter_parameter_mpt
  use interface_violation_midpoint_MoC_CF_peos_spin
  use interface_violation_gridpoint_MoC_CF_peos_spin
  use interface_violation_midpoint_HaC_CF_peos_spin
  use interface_violation_gridpoint_HaC_CF_peos_spin
  use interface_IO_output_3D_general
  use interface_IO_output_2D_general
  use interface_IO_output_1D_general
  use make_array_2d
  use make_array_3d
  use make_array_4d
  implicit none
  real(long), pointer :: pot(:,:,:), HaC_vio(:,:,:), MoC_vio(:,:,:,:)
  real(long) :: xa,xb,fa,fb,rm_eps,delta,f1,f2, iter_eps
  integer :: total_iteration, iter_count, i, icycle, Nmax=15
  integer :: impt
  character(30) :: char1, char2, char3, char4, char5

  rm_eps = 1.0d-06
  total_iteration = 0
  char3 = 'main_bnsphys_all_mpt.txt'
!
  call open_directory_mpt(0)
  iter_eps = 1.0d0
  call iter_lecc_spin_BNS_CF_mpt(iter_count,1,iter_eps)     !  1 for freezing hydro in first 5 iterations
  total_iteration = total_iteration + iter_count
  write(6,'(a21,i5)') '--- total iteration =', total_iteration
  call calc_physical_quantities_spin_BNS_CF_mpt   
  call write_last_physq_BNS_mpt
  call printout_physq_BNS_all_mpt(char3)
!  call save_solution_BNS_mpt(0)
  call copy_def_matter_parameter_from_mpt(1)
  xa = emdc
  call copy_def_quantities_BNS_from_mpt(1)
  fa = (restmass - restmass_sph)/restmass_sph
  write(6,*)  "xa,fa=", xa, fa,  "           restmass, restmass_sph = ", restmass, restmass_sph
  call chdir('../')
!
  call open_directory_mpt(1)
!  xb = xa - 0.25d0*xa
  xb = xa - 0.01d0*xa
  do impt = 1, 2
    call copy_def_matter_parameter_from_mpt(impt)
    emdc = xb
    call copy_def_matter_parameter_to_mpt(impt)
  end do
  write(6,'(a14,1p,2e20.12)')  "emdc COCP1,2: ", def_matter_param_real_(2,1), def_matter_param_real_(2,2)
  iter_eps = 1.0d0
  call iter_lecc_spin_BNS_CF_mpt(iter_count,2,iter_eps)     ! 2 to compute hydro from first iteration
  total_iteration = total_iteration + iter_count
  write(6,'(a21,i5)') '--- total iteration =', total_iteration
  call calc_physical_quantities_spin_BNS_CF_mpt
  call write_last_physq_BNS_mpt
  call printout_physq_BNS_all_mpt(char3)
!  call save_solution_BNS_mpt(1)
  call copy_def_quantities_BNS_from_mpt(1)
  fb = (restmass - restmass_sph)/restmass_sph
  write(6,*)  "xb,fb=", xb, fb,  "           restmass, restmass_sph = ", restmass, restmass_sph
  call chdir('../')
!
  write(6,*) '--------- Initial values for emdc and f(emdc)=(restmass - restmass_sph)/restmass_sph ---------------------------'
  write(6,'(a30,i2,1p,2e20.12)')  '(cycle,emdc,(M0-M0sph)/M0sph)=', 0,xa,fa
  write(6,'(a30,i2,1p,2e20.12)')  '(cycle,emdc,(M0-M0sph)/M0sph)=', 1,xb,fb
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
    if ( dabs(delta) < rm_eps ) then
      write(6,*) 'Convergence...'
      call copy_def_matter_parameter_from_mpt(1)
      write(6,*) '------------------------------------------------------------------'
      write(6,*) 'Central Emden = ',emdc

!      write(6,*) "Calculating constraint violations at gridpoints of COCP-1..."
!      call copy_from_mpatch_all_BNS_CF(1)
!      call copy_def_metric_and_matter_from_mpt(1)
!      call copy_def_matter_spin_from_mpt(1)
!      call alloc_array4d(MoC_vio,0,nrg,0,ntg,0,npg,1,3)
!      call alloc_array3d(HaC_vio,0,nrg,0,ntg,0,npg)
!      call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
!      call excurve_CF('ns')             !   3rd order midpoint from ir0=-2,...
!      call excurve_CF_gridpoint         !   4th order from ir0=-2,...

!      MoC_vio(0:nrg,0:ntg,0:npg,1:3) = 0.0d0
!      call violation_gridpoint_MoC_CF_peos_spin(MoC_vio,'ns')
!      pot = 0.0d0
!      pot(0:nrg,0:ntg,0:npg) = MoC_vio(0:nrg,0:ntg,0:npg,2)
!      char1 = 'MoC_by_3D_mpt1.txt'
!      call IO_output_3D_general(char1,'g','g',pot)
!      char1 = 'MoC_by_xy_mpt1.txt'
!      call IO_output_2D_general(char1,'g','g',pot,'xy')
!      char1 = 'MoC_by_phi000_mpt1.txt'
!      call IO_output_1D_general(char1,'g','g',pot,-1,ntg/2,0)
!      char1 = 'MoC_by_phi180_mpt1.txt'
!      call IO_output_1D_general(char1,'g','g',pot,-1,ntg/2,npg/2)

!      HaC_vio = 0.0d0
!      call violation_gridpoint_HaC_CF_peos_spin(HaC_vio,'ns')
!      char1 = 'HaC_3D_mpt1.txt'
!      call IO_output_3D_general(char1,'g','g',HaC_vio)
!      char1 = 'HaC_xy_mpt1.txt'
!      call IO_output_2D_general(char1,'g','g',HaC_vio,'xy')
!      char1 = 'HaC_phi000_mpt1.txt'
!      call IO_output_1D_general(char1,'g','g',HaC_vio,-1,ntg/2,0)
!      char1 = 'HaC_phi180_mpt1.txt'
!      call IO_output_1D_general(char1,'g','g',HaC_vio,-1,ntg/2,npg/2)

!      deallocate(MoC_vio)
!      deallocate(HaC_vio)
!      deallocate(pot)
!      call printout_physq_console_BBH
      return
    end if
!    xa = xa - 0.1d0*delta
    if (i.lt.10) then
      xa = xa - i*delta/10.0d0
    else
      xa = xa - delta
    end if 
    if (i.lt.7) then
      iter_eps = dmin1(10.0d0**(-i), dabs(delta))
    else
      iter_eps = rm_eps
    end if

    call open_directory_mpt(icycle)
    do impt = 1, 2
      call copy_def_matter_parameter_from_mpt(impt)
      emdc = xa
      call copy_def_matter_parameter_to_mpt(impt)
    end do
    write(6,'(a14,1p,2e20.12)')  "emdc COCP1,2: ", def_matter_param_real_(2,1), def_matter_param_real_(2,2)
    write(6,*) '------------------------------------------------------------------'
!    write (6,'(a7,i2,a20,e20.12)') 'cycle =', icycle, '        New emdc  = ',emdc

    call iter_lecc_spin_BNS_CF_mpt(iter_count,icycle,iter_eps)
    total_iteration = total_iteration + iter_count
    write(6,'(a21,i5)') '--- total iteration =', total_iteration
    call calc_physical_quantities_spin_BNS_CF_mpt
    call write_last_physq_BNS_mpt
    call printout_physq_BNS_all_mpt(char3)
!    call save_solution_BNS_mpt(icycle)
    call copy_def_quantities_BNS_from_mpt(1)
    fa = (restmass - restmass_sph)/restmass_sph
    write(6,'(a30,i2,1p,2e20.12)')   '(cycle,emdc,(M0-M0sph)/M0sph)=',icycle,xa,fa
    call chdir('../')
!
  end do
  write(6,*) 'No convergence after ',Nmax,' iterations.'

end subroutine calc_lecc_spin_BNS_CF_mpt
