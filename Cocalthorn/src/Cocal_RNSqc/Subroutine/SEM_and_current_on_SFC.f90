subroutine SEM_and_current_on_SFC
  use phys_constant, only :  long
  implicit none
  real(long) :: tic1, tic2
!
  call interpo_gr2fl_metric_WL
  call calc_h_udf_axisym
  call SEM_tensor
  call calc_Aphi_max
!!  call calc_MHDpar_charge
  call current_jt_MHD('ns')
  call current_MHD_on_SFC
  call current_sfc2sph_MHD
!
!!call cpu_time(tic1)
!!call cpu_time(tic2)
!!write(6,*)'time current_jt_MHD  ', tic2-tic1 ; tic1=tic2
!  call current_MHD
!  call current_sfc2sph_MHD
!  call current_modify_MHD('ns')
!  call current_sph2sfc_MHD
end subroutine SEM_and_current_on_SFC
