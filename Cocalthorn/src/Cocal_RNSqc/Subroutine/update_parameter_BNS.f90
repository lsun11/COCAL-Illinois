subroutine update_parameter_BNS(conv_den)
  use grid_parameter, only : EQ_point, nrf_deform, nrf, NS_shape, &
                        &    sw_sepa
  use def_matter_parameter, only : ROT_LAW
  use interface_update_parameter_BNS_peos
  use interface_update_parameter_BNS_peos_lecc
  use interface_update_parameter_BNS_peos_irrot
  use interface_update_parameter_BNS_peos_lecc_irrot
  use interface_update_parameter_BNS_peos_spin
  use interface_update_parameter_BNS_peos_lecc_spin
!  use interface_update_parameter_axisym_peos
!  use interface_update_parameter_triaxial_peos
!  use interface_update_parameter_spherical_peos
  use interface_update_parameter_axisym_peos_drot
  implicit none
  real(8), intent(in) :: conv_den
  if (ROT_LAW.eq.'DR') then 
    call update_parameter_axisym_peos_drot(conv_den)
  else 
!    call update_parameter_BNS_peos(conv_den)

    if (NS_shape.eq.'CO')   then
      if (EQ_point.eq.'LE')  then
        call update_parameter_BNS_peos_lecc(conv_den)
      else 
        call update_parameter_BNS_peos(conv_den)
      end if
    end if

    if (NS_shape.eq.'IR')   then
      if (EQ_point.eq.'LE')  then
        call update_parameter_BNS_peos_lecc_irrot(conv_den)
      else 
        call update_parameter_BNS_peos_irrot(conv_den)
      end if
    end if

    if (NS_shape.eq.'SP')   then
      if (EQ_point.eq.'LE')  then
        call update_parameter_BNS_peos_lecc_spin(conv_den)
      else 
        call update_parameter_BNS_peos_spin(conv_den)
      end if
    end if

    if(NS_shape.ne.'CO' .and. NS_shape.ne.'IR' .and. NS_shape.ne.'SP')  then
      write(6,*) "NS_shape is neither CO nor IR nor SP. Exiting..."
      stop
    end if
!    if (EQ_point.eq.'XZ'.and.nrf.ne.nrf_deform) &
!    &                     call update_parameter_axisym_peos(conv_den)
!    if (EQ_point.eq.'XY') call update_parameter_triaxial_peos(conv_den)
!    if (EQ_point.eq.'XZ'.and.nrf.eq.nrf_deform) &
!    &                     call update_parameter_spherical_peos(conv_den)
  end if
end subroutine update_parameter_BNS
