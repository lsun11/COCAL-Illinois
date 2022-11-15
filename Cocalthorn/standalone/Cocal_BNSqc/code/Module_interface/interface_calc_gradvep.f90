module interface_calc_gradvep
  implicit none
  interface 
    subroutine calc_gradvep(potf,potxf,potyf,potzf)
      real(8), pointer ::  potf(:,:,:), potxf(:,:,:), potyf(:,:,:), potzf(:,:,:)
    end subroutine calc_gradvep
  end interface
end module interface_calc_gradvep
