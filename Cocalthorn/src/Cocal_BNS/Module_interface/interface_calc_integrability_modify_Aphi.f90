module interface_calc_integrability_modify_Aphi
  implicit none
  interface 
    subroutine calc_integrability_modify_Aphi(fncx,fncy,supp)
      real(8), pointer :: fncx(:,:,:), fncy(:,:,:)
      character(len=2), intent(in) :: supp
    end subroutine calc_integrability_modify_Aphi
  end interface
end module interface_calc_integrability_modify_Aphi
