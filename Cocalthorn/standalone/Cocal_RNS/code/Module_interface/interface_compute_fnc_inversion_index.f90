module interface_compute_fnc_inversion_index
  implicit none
  interface
    subroutine compute_fnc_inversion_index(pot1,pot2,index)
      real(8), pointer :: pot1(:,:,:), pot2(:,:,:)
      real(8) :: index
    end subroutine compute_fnc_inversion_index
  end interface
end module interface_compute_fnc_inversion_index
