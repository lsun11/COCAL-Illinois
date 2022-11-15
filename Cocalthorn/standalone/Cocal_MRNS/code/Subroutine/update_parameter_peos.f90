subroutine update_parameter_peos(conv_den)
  use grid_parameter, only : EQ_point, nrf_deform, nrf
  use def_matter_parameter, only : ROT_LAW
  use interface_update_parameter_axisym_peos
  use interface_update_parameter_triaxial_peos
  use interface_update_parameter_spherical_peos
  use interface_update_parameter_axisym_peos_drot
  implicit none
  real(8), intent(in) :: conv_den
!
  if (nrf.ne.nrf_deform) then 
    if (ROT_LAW.eq.'DR'.or.ROT_LAW.eq.'OJ') then 
      call update_parameter_axisym_peos_drot(conv_den)
    else if (EQ_point.eq.'XZ') then
      call update_parameter_axisym_peos(conv_den)
    else if (EQ_point.eq.'XY') then 
      call update_parameter_triaxial_peos(conv_den)
    end if
  else
    call update_parameter_spherical_peos(conv_den)
  end if
end subroutine update_parameter_peos
