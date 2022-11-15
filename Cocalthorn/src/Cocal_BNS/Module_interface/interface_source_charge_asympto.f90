module interface_source_charge_asympto
  implicit none
  interface 
    subroutine source_charge_asympto(sousf,irg)
      real(8), pointer :: sousf(:,:)
      integer          :: irg
    end subroutine source_charge_asympto
  end interface
end module interface_source_charge_asympto
