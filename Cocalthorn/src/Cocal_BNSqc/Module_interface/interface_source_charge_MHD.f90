module interface_source_charge_MHD
  implicit none
  interface 
    subroutine source_charge_MHD(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine source_charge_MHD
  end interface
end module interface_source_charge_MHD
