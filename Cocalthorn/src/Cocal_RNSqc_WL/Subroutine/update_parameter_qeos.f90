subroutine update_parameter_qeos(conv_den)
  use grid_parameter, only : EQ_point, nrf_deform, nrf
  use def_matter_parameter, only : ROT_LAW
  use interface_update_parameter_axisym_qeos
  use interface_update_parameter_triaxial_qeos
  use interface_update_parameter_spherical_qeos
  use interface_update_parameter_axisym_qeos_drot
  implicit none
  real(8), intent(in) :: conv_den
  if (ROT_LAW.eq.'DR') then 
    call update_parameter_axisym_qeos_drot(conv_den)
  else 
    if (EQ_point.eq.'XZ'.and.nrf.ne.nrf_deform) &
    &                     call update_parameter_axisym_qeos(conv_den)
    if (EQ_point.eq.'XY') call update_parameter_triaxial_qeos(conv_den)
    if (EQ_point.eq.'XZ'.and.nrf.eq.nrf_deform) &
    &                     call update_parameter_spherical_qeos(conv_den)
  end if
end subroutine update_parameter_qeos
