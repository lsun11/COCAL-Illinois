subroutine update_parameter(conv_den)
  use grid_parameter, only : EQ_point
  use interface_update_parameter_axisym
  use interface_update_parameter_triaxial
  implicit none
  real(8), intent(in) :: conv_den
  if (EQ_point.eq.'XZ') call update_parameter_axisym(conv_den)
  if (EQ_point.eq.'XY') call update_parameter_triaxial(conv_den)
end subroutine update_parameter
