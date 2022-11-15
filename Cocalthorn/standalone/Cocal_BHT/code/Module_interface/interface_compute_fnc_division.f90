module interface_compute_fnc_division
  implicit none
  interface 
    subroutine compute_fnc_division(pot1,pot2,pot3)
      real(8), pointer :: pot1(:,:,:), pot2(:,:,:), pot3(:,:,:)
    end subroutine compute_fnc_division
  end interface
end module interface_compute_fnc_division
