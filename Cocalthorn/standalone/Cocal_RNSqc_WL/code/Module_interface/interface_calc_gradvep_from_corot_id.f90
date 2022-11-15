module interface_calc_gradvep_from_corot_id
  implicit none
  interface 
    subroutine calc_gradvep_from_corot_id(potf,potxf,potyf,potzf)
      real(8), pointer ::  potf(:,:,:), potxf(:,:,:), potyf(:,:,:), potzf(:,:,:)
    end subroutine calc_gradvep_from_corot_id
  end interface
end module interface_calc_gradvep_from_corot_id
