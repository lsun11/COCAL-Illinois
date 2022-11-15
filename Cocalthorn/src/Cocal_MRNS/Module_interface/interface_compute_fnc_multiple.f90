module interface_compute_fnc_multiple
  implicit none
  interface
    subroutine compute_fnc_multiple(pot1,pot2,pot3)
      real(8), pointer :: pot1(:,:,:), pot2(:,:,:), pot3(:,:,:)
    end subroutine compute_fnc_multiple
  end interface
end module interface_compute_fnc_multiple
