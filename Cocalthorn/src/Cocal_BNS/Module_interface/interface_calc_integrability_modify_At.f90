module interface_calc_integrability_modify_At
  implicit none
  interface 
    subroutine calc_integrability_modify_At(fnc,supp)
      real(8), pointer :: fnc(:,:,:)
      character(len=2), intent(in) :: supp
    end subroutine calc_integrability_modify_At
  end interface
end module interface_calc_integrability_modify_At
